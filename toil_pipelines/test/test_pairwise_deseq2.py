import os
from subprocess import check_call


def test_workflow(tmpdir):
    workdir = str(tmpdir)

    # Specify path to sample dataframe which is 20 genes by 20 samples
    cwd = os.path.dirname(os.path.abspath(__file__))
    df_path = os.path.join(cwd, 'sample-20-x-20.tsv')
    assert os.path.exists(df_path)

    # Create sample pairing information
    samples = open(df_path, 'r').readline().split()
    assert len(samples) == 20
    pairing_path = os.path.join(workdir, 'pairing.tsv')
    with open(pairing_path, 'w') as f:
        for pair in zip(samples, ['1'] * 10 + ['2'] * 10):
            f.write('\t'.join(pair) + '\n')

    # Create manifest
    manifest_path = os.path.join(workdir, 'manifest-pairwise-deseq2.tsv')
    with open(manifest_path, 'w') as f:
        f.write('test\tfile://{}\tfile://{}'.format(df_path, pairing_path))

    # Call pipeline
    pipeline_path = os.path.join(os.path.dirname(cwd), 'pairwise_deseq2.py')
    check_call(['python', pipeline_path,
                'run',
                '--manifest', manifest_path,
                '--output-dir', workdir,
                '--initial-size', '1G',
                '--max-partition', '5',
                '--retryCount', '2',
                '--workDir', workdir,
                os.path.join(workdir, 'jobStore')])

    # Confirm output exists
    assert os.path.exists(os.path.join(workdir, 'test.tsv'))


