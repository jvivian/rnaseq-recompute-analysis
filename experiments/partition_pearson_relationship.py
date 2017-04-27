"""
Attempt to ascertain the tradeoff between:
    - The size of the partition
    - The # of times it is run

Method:
For a given partition size
    Run Pairwise DESeq2
    Save all Results
    Save jobstore to query stats
"""
import os
import shutil
from subprocess import check_call, CalledProcessError

from experiments.AbstractExperiment import AbstractExperiment


class PartitionPearsonRelationship(AbstractExperiment):
    def __init__(self, root_dir):
        super(PartitionPearsonRelationship, self).__init__(root_dir)
        self.experiment_dir = os.path.join(root_dir, 'experiments/partition-pearson-relationship')
        self.manifest_dir = os.path.join(self.experiment_dir, 'manifests')
        self.log_dir = os.path.join(self.experiment_dir, 'logs')
        self.results_dir = os.path.join(self.experiment_dir, 'results')
        self.work_dir = os.path.join(self.experiment_dir, 'workDir')
        self.jobStores = os.path.join(self.experiment_dir, 'jobStores')

        # Tunable parameters
        self.partitions = ['8']  # , 2, 4, 8, 16, 32, 64, 128, 273, 600, 1092]
        self.tissue = 'breast'
        self.expression_path = 'file://' + [x for x in self.protein_coding_paths if 'breast' in x][0]
        self.group_path = 's3://jvivian-recompute-analysis/combat-corrected-dataframes/breast/tcga.tsv'

    def setup(self):
        dirtree = [os.path.join(self.results_dir, p) for p in self.partitions] + \
                  [self.manifest_dir, self.work_dir, self.jobStores, self.log_dir]
        self.create_directories(dirtree)

    def run_experiment(self):
        cmd = ['time', 'pairwise-deseq2', 'run', '--retryCount=2']
        for p in self.partitions:

            manifest_path = os.path.join(self.manifest_dir, 'manifest-{}.tsv'.format(p))
            with open(manifest_path, 'w') as f:
                f.write('{}\t{}\t{}\n'.format('{}-{}'.format(self.tissue, p), self.expression_path, self.group_path))

            output_dir = os.path.join(self.results_dir, p)
            with open(os.path.join(self.log_dir, 'log-{}.txt'.format(p)), 'w') as f:
                check_call(cmd + ['--manifest={}'.format(manifest_path),
                                  '--max-partition=' + str(p),
                                  os.path.join(self.jobStores, str(p)),
                                  '--workDir='.format(os.path.join(self.work_dir)),
                                  '--output-dir=' + output_dir,
                                  '--disableCaching',
                                  '--cleanWorkDir=never',
                                  '--clean=never',
                                  '--stats'],
                           stdout=f, stderr=f)

            # Collect Results
            toil_workdir = os.listdir(self.work_dir)[0]
            i = 0
            for root, d, files in os.walk(toil_workdir):
                if 'results.tsv' in files:
                    shutil.move(os.path.join(root, 'results.tsv'),
                                os.path.join(self.results_dir, p, str(i) + '.tsv'))
                    i += 1
            shutil.rmtree(toil_workdir)

    def teardown(self):
        pass
