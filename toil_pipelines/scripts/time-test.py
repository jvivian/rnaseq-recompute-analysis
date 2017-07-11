import os
from subprocess import Popen, PIPE

# This script is to produce data for a single plot, so paths are hardcoded
df_path = '/mnt/rna-seq-analysis/data/xena/tcga_gtex_counts_protein_coding.tsv'
pairing_path = '/home/ubuntu/rnaseq-recompute-analysis/toil_pipelines/scripts/pairing.tsv'
cwd = os.path.dirname(os.path.abspath(__file__))
pipeline_path = os.path.join(os.path.dirname(cwd), 'pairwise_deseq2.py')
manifest_path = os.path.join(cwd, 'manifest-pairwise-deseq2.tsv')
log_path = os.path.join(cwd, 'log.txt')

# Delete log if exists as we'll be appending
if os.path.exists(log_path):
    os.remove(log_path)

# For every log2 partition run DESeq2 and save output
for max_partition in [32, 64, 128, 256, 512, 1024]:

    # Create manifest
    with open(manifest_path, 'w') as f:
        f.write('test{}\tfile://{}\tfile://{}'.format(max_partition, df_path, pairing_path))

    print 'Running workflow with a maximum partition of: {}'.format(max_partition)
    # Call workflow
    p = Popen(['time', 'python', pipeline_path,
               'run',
               '--manifest', manifest_path,
               '--output-dir', cwd,
               '--initial-size', '1G',
               '--max-partition', str(max_partition),
               '--retryCount', '2',
               '--workDir', '/mnt/',
               '/mnt/jobStore'], stderr=PIPE, stdout=PIPE)

    out, err = p.communicate()

    if p.returncode != 0:
        print 'Something went FUBAR with run: {}'.format(max_partition)
        print out
        print err
    else:
        with open(log_path, 'a') as f:
            f.write('{}\n\n'.format(err))
