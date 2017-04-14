from __future__ import print_function

import argparse
import os
import pickle
import sys
import textwrap
from collections import defaultdict
from subprocess import check_call
from urlparse import urlparse
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.common import Toil
from toil.job import Job, PromisedRequirement
from toil_lib import require, UserError, flatten
from toil_lib.files import copy_files
from toil_lib.urls import download_url, s3am_upload, download_url_job


schemes = ('http', 'file', 's3', 'ftp', 'gnos')


def root(job, samples, gene_map, output_dir, disk):
    """
    Downloads input file and stages each sample

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param dict samples: key=UUID, value=[expression_url, group_url]
    :param str gene_map: URL to gene_map used to map gene IDs to gene names
    :param str|int disk:
    :param str output_dir: Full path/ S3 URL to output directory
    """
    # Download gene_map
    gene_map_id = job.addChildJobFn(download_url_job, gene_map).rv()

    # Staging is IO intensive, and should run one job per machine if possible
    cores = int(cpu_count())

    # For every sample -> staging
    for sample in samples.iteritems():
        job.addFollowOnJobFn(staging, sample, gene_map_id, output_dir, disk=disk, cores=cores).rv()


def staging(job, sample, gene_map_id, output_dir):
    """
    Derives vectors for a sample and spawns jobs for individual DESeq2 runs

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param tuple(str, [str, str]) sample: tuple containing UUID and URLs for the expression dataframe and group info
    :param str gene_map_id: FileStoreID for the gene_map downloaded in the root job
    :param str output_dir: Full path/ S3 URL to output directory
    """
    work_dir = job.fileStore.getLocalTempDir()
    uuid, df_url, group_url = flatten(sample)

    # Download dataframe and group info
    group_path = download_url(group_url, work_dir=work_dir)
    df_path = download_url(df_url, work_dir=work_dir)
    df_id = job.fileStore.writeGlobalFile(df_path)

    # Parse group information
    groups = defaultdict(list)
    with open(group_path, 'r') as f:
        for line in f:
            if not line.isspace():
                sample, group = line.strip().split('\t')
                groups[group].append(sample)
    assert len(groups.keys()) == 2 and '1' in groups.keys() and '2' in groups.keys(), \
        'Group requirements not met. Make sure only a 1 or a 2 is used: {}'.format(groups.keys())

    # Sort so group 1 is smaller
    group_1, group_2 = sorted([groups['1'], groups['2']], key=lambda x: len(x))

    log(job, 'Writing out vectors, one per sample in group 2. g1: {}  g2: {}'.format(len(group_1), len(group_2)))
    for sample in group_2:
        with open(os.path.join(work_dir, sample), 'w') as f:
            vector = [x.replace('-', '.') for x in group_1 + [sample]]  # Precaution for R
            f.write('\n'.join(vector) + '\n')

    # Save vectors to jobStore
    vector_ids = [job.fileStore.writeGlobalFile(os.path.join(work_dir, x)) for x in group_2]

    log(job, 'Creating one job per vector: {}'.format(len(vector_ids)))
    disk = int(df_id.size + 1e9)
    results = [job.addChildJobFn(run_deseq2, df_id, vector_id, disk=disk).rv() for vector_id in vector_ids]

    # Follow-on --> combine_results, return
    disk = PromisedRequirement(lambda y: sum([x.size for x in y]) + int(1e9), results)
    # disk = int(len(group_2) * 1e9)
    job.addFollowOnJobFn(combine_results, results, gene_map_id, uuid, output_dir, disk=disk).rv()


def run_deseq2(job, df_id, vector_id):
    """
    Runs DEseq2 given a vector and expression dataframe

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param str df_id: FileStoreID of expression dataframe
    :param str vector_id: FileStoreID of vector
    :return: FileStoreID of the results
    :rtype: str
    """
    # Read in inputs
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.readGlobalFile(df_id, os.path.join(work_dir, 'expression.tsv'))
    job.fileStore.readGlobalFile(vector_id, os.path.join(work_dir, 'vector'))

    # Write out Rscript
    r_path = os.path.join(work_dir, 'deseq.R')
    with open(r_path, 'wb') as f:
        f.write(deseq2_script())

    # Run DESeq2 Docker
    tool = 'jvivian/deseq2:1.14.1'
    base_command = ['docker', 'run',
                    '--rm',
                    '--log-driver=none',
                    '-v', '{}:/data'.format(work_dir)]

    deseq2_command = base_command + [tool, '/data/deseq.R']
    check_call(deseq2_command)

    # Fix permissions on output by reusing tool container
    stat = os.stat(work_dir)
    check_call(base_command + ['--entrypoint=chown', tool, '-R', '{}:{}'.format(stat.st_uid, stat.st_gid), '/data'])

    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'results.tsv'))


def combine_results(job, result_ids, gene_map_id, uuid, output_dir):
    """
    Combines individual DESeq2 results into one meta-results TSV

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list[str,] result_ids: FileStoreIDs of all the results for this sample
    :param str gene_map_id: FileStoreID for the gene_map downloaded in the root job
    :param str uuid: The UUID of the sample
    :param str output_dir: Full path/ S3 URL to output directory
    """
    log(job, 'Combining results for: ' + uuid)
    work_dir = job.fileStore.getLocalTempDir()
    gene_map_path = job.fileStore.readGlobalFile(gene_map_id, os.path.join(work_dir, 'gene_map.pickle'))
    gene_map = pickle.load(open(gene_map_path, 'rb'))

    # Read in results
    results = []
    for i, result_id in enumerate(result_ids):
        if result_id is not None:
            results.append(job.fileStore.readGlobalFile(result_id, os.path.join(work_dir, '{}.tsv'.format(i))))

    # Collect results for each gene
    ranked = pd.DataFrame()
    pvals = defaultdict(list)
    fc = defaultdict(list)
    for result in results:
        with open(result, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line:
                    gene, mean, l2fc, lfcse, stat, pval, padj = line
                    if padj != 'NA':
                        pvals[gene].append(float(padj))
                        fc[gene].append(float(l2fc))

    # Create columns in the combined dataframe
    genes = pvals.keys()
    ranked['num_samples'] = [len(pvals[x]) for x in genes]
    ranked['pval_counts'] = [sum([1 for y in pvals[x] if y < 0.01]) for x in genes]
    ranked['pval_percentage'] = np.array([float(x) for x in ranked['pval_counts']]) / np.array(ranked['num_samples'])
    ranked['pval'] = [np.median(pvals[x]) for x in genes]
    ranked['pval_std'] = [round(np.std(pvals[x]), 4) for x in genes]
    ranked['fc'] = [round(np.median(fc[x]), 4) for x in genes]
    ranked['fc_std'] = [round(np.std(fc[x]), 4) for x in genes]

    # Map gene ids to gene names and add as a column
    gene_names = [gene_map[x] if x in gene_map.keys() else x for x in genes]
    ranked['gene_id'] = genes
    ranked.index = gene_names
    ranked.sort_values('pval_counts', inplace=True, ascending=False)

    # Write out dataframe
    results_path = os.path.join(work_dir, '{}.tsv'.format(uuid))
    ranked.to_csv(results_path, sep='\t')

    # Upload or move to directory
    if urlparse(output_dir).scheme == 's3':
        log(job, 'Uploading {} to S3: {}'.format(uuid, output_dir))
        s3am_upload(fpath=results_path, s3_dir=output_dir)
    else:
        log(job, 'Moving {} to output dir: {}'.format(uuid, output_dir))
        mkdir_p(output_dir)
        copy_files(file_paths=[results_path], output_dir=output_dir)


# Pipeline specific functions
def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each analysis to run
        #   There are 3 tab-separated columns: UUID, dataframe URL, and group info URL.
        #   Valid URL schemes: {scheme}
        #
        #   UUID        Unique identifier for the sample
        #
        #   dataframe   Tab-separated gene by sample matrix
        #
        #   group       A two column TSV assocating sample to group number. First column is a sample
        #               name corresponding to a sample in the dataframe, and the second column is
        #               either a 1 or a 2 indicating which group the sample is associated with.
        #
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/expression.tsv  file:///path/to/group_info.tsv
        #   UUID_2  http://sample-depot.com/expression.tsv   file:///path/to/group_info.tsv
        #   UUID_3  s3://my-bucket-name/dir/expression.tsv  s3://my-bucket-name/dir/group_info.tsv
        #
        #   Place your samples below, one per line.
        """.format(scheme=[x + '://' for x in schemes])[1:])


def parse_samples(path_to_manifest):
    """
    Parses samples as specified in the manifest

    :param str path_to_manifest: Path to manifest
    :return: Sample information from manifest
    :rtype: dict
    """
    samples = {}
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split('\t')
                require(len(sample) == 3, 'Bad manifest format! '
                                          'Expected 3 tab separated columns, got: {}'.format(sample))
                uuid, exp_url, group_url = sample
                for url in [exp_url, group_url]:
                    require(urlparse(url).scheme in schemes, 'Samples must start with one of the '
                                                             'approved schemes: {}'.format(schemes))
                samples[uuid] = [exp_url, group_url]
    return samples


def deseq2_script():
    return textwrap.dedent("""
        suppressMessages(library('DESeq2'))
        suppressMessages(library('data.table'))

        # Argument parsing
        df_path <- '/data/expression.tsv'
        vector_path <- '/data/vector'
        results_dir <- '/data'
        sample_name <- 'results.tsv'

        # Read in tables / patients
        n <- read.table(df_path, sep='\\t', header=1, row.names=1)
        vector <- read.table(vector_path)$V1
        sub <- n[, colnames(n)%in%vector]
        setcolorder(sub, as.character(vector))

        # Create matrix vectors
        disease_vector <- c(rep('A', length(vector)-1), 'B')

        # DESeq2 preprocessing
        # Rounding the countData since DESeQ2 only accepts integer counts
        countData <- round(sub)
        colData <- data.frame(disease=disease_vector, row.names=colnames(countData))
        y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ disease)

        # Run DESeq2
        y <- DESeq(y)
        res <- results(y)
        summary(res)

        # Write out table
        resOrdered <- res[order(res$padj),]
        res_path <- paste(results_dir, sample_name, sep='/')
        write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\\t', quote=FALSE)
        """)


def log(job, string):
    job.fileStore.logToMaster('\n\n\n{}\n\n'.format(string))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Pairwise DESeq2 Pipeline

    RNA-seq expression data is run pairwise and combined using DESeq2

    General usage:
    1. Type "pairwise-deseq2 generate" to create an editable manifest and config in the current working directory.
    2. Fill in the manifest with information pertaining to your samples.
    4. Type "pairwise-deseq2 run [jobStore]" to execute the pipeline.

    Please read the README.md located in the source directory or at:

    Structure of the Pairwise DESeq2 Pipeline

        n        m
    0 -----> 1 ------> 2
    |        |
    |        |
    4        3

    0 = Root job. Runs n jobs per samples submitted
    1 = Sample staging. Runs m jobs per samples in the larger group.
    2 = Runs DESeq2 comparing the samller group to a single sample from the larger group
    3 = Combines results into a dataframe sorting genes by p-value counts. Moves results to output-dir
    =======================================
    Dependencies
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL, needed for upload to S3)
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    # Generate subparser
    subparsers.add_parser('generate', help='Generates an editable manifest in the current working directory.')
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the RNA-seq pipeline')
    parser_run.add_argument('--manifest', default='manifest-pairwise-deseq2.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')

    parser_run.add_argument('--output-dir', required=True, help='Output location for pipeline. Either provide '
                                                                'a full path to a directory or an s3:// URL '
                                                                'where the samples will be uploaded.')

    parser_run.add_argument('--gene_map', help='A python pickled dictionary mapping gene ids to gene names',
                            default='https://raw.githubusercontent.com/jvivian/'
                                    'rnaseq-recompute-analysis/master/utils/gene_map.pickle')
    parser_run.add_argument('--initial-size', default='10G', help='Initial disk allotted for each sample. Change'
                                                                  'this value if your average sample is larger'
                                                                  'than 10G.')
    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()
    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-pairwise-deseq2.tsv'), generate_manifest)
    # Pipeline execution
    elif args.command == 'run':
        samples = parse_samples(path_to_manifest=args.manifest)
        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        with Toil(args) as toil:
            if args.restart:
                toil.restart()
            else:
                toil.start(Job.wrapJobFn(root, samples, args.gene_map, args.output_dir, args.initial_size))


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
