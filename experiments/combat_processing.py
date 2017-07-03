import logging
import os
import shutil
import textwrap
from subprocess import check_call

import pandas as pd

from experiments.abstractExperiment import AbstractExperiment
from utils import quantile_normalize
from utils import write_script

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class CombatProcessing(AbstractExperiment):

    def __init__(self, root_dir, cores, parametric, tissue_cov, s3_dir):
        super(CombatProcessing, self).__init__(root_dir)
        self.cores = cores
        self.experiment_dir = os.path.join(self.experiments_dir, '')
        self.combat_tsv = 'expression-CB'

        # Parameters for ComBat - edit as needed
        self.parametric = parametric  # Boolean
        self.tissue_cov = tissue_cov  # Boolean
        self.s3_dir = s3_dir  # String or None
        if self.s3_dir:
            assert self.s3_dir.startswith('s3://'), 'S3 format incorrect. Must start with "s3://"'

        # Change output name contingent upon parameters
        self.combat_tsv += 'P-' if self.parametric else 'NP-'
        self.combat_tsv += 'T' if self.tissue_cov else ''
        self.combat_tsv = os.path.join(self.experiment_dir, self.combat_tsv + '.tsv')

    def setup(self):
        pass

    def run_experiment(self):
        log.info('Creating master dataframe')
        df = pd.concat([pd.read_csv(t, sep='\t', index_col=0)for t in self.protein_coding_paths], axis=1)

        log.info('Dropping duplicate columns by name')
        df = df.T.groupbyu(level=0).first().T

        # Quantile Normalization
        g = [x for x in df.columns if x.startswith('GTEX')]
        t = [x for x in df.columns if x.startswith('TCGA') and (x.endswith('11') or x.endswith('01'))]
        log.info('Quantile Normalizing GTEx')
        gtex = quantile_normalize(df[g])
        log.info('Quantile Normalizing TCGA')
        tcga = quantile_normalize(df[t])

        log.info('Combining separately quantile-normalized dataframes')
        df = pd.concat([gtex, tcga], axis=1)

        log.info('Removing low variance genes')
        df = df[df.std(axis=1) > 1]

        log.info('Writing out dataframe')
        df.to_csv(self.combat_tsv, sep='\t')

        log.info('Creating Sample to tissue mapping')
        tmap = {}
        for t in self.protein_coding_paths:
            with open(t, 'r') as f_in:
                samples = f_in.readline().strip().split('\t')
            for sample in samples:
                tmap[sample] = t

        log.info('Writing out batch.tsv')
        with open(self.experiment_dir, 'batch.txt') as f:
            f.write('id\tBatch\tCovariate\n') if self.tissue_cov else f.write('id\tBatch\n')
            for sample in df.columns:
                batch = 'gtex' if sample.startswith('GTEX') else 'tcga'
                if self.tissue_cov:
                    f.write('{}\t{}\t{}\n'.format(sample, batch, tmap[sample]))
                else:
                    f.write('{}\t{}\n'.format(sample, batch))

        log.info('Running ComBat')
        shutil.copy(self.combat_script, self.experiment_dir)
        write_script(self.run_combat_script(), self.experiment_dir, name='run_combat.r')
        check_call(['Rscript', os.path.join(self.experiment_dir, 'combat.r')])

        log.info('Reading in ComBat corrected dataframe')
        output = os.path.join(self.experiment_dir, 'Adjusted_' + os.path.basename(self.combat_tsv) + '_.xls')
        df = pd.read_csv(output, sep='\t', index_col=0)

        log.info('Dropping ProbeID')
        df = df.drop('ProbeID', axis=1)

        log.info('Removing negative values from expression array')
        df[df < 0] = 0

        log.info('Adding back in genes')
        genes = []
        with open(self.combat_tsv, 'r') as f:
            f.readline()
            for line in f:
                genes.append(line.split('\t')[0])
        df.index = genes

        log.info('Subsetting Dataframe by tissue')
        base = ['s3am', 'upload', '--exists=overwrite']
        for t in self.protein_coding_paths:
            tissue = os.path.basename(os.path.dirname(t))
            log.info(os.path.basename(t))
            with open(t, 'r') as f:
                samples = f.readline().strip().split('\t')

            samples = [x for x in samples if x.startswith('GTEX') or (x.endswith('11') or x.endswith('01'))]
            temp = df[samples]

            f = os.path.join(self.tissue_pair_dir, tissue, os.path.basename(self.combat_tsv))
            temp.to_csv(f, sep='\t')

            if self.s3_dir:
                log.info('Uploading to S3: ' + self.s3_dir)
                check_call(base + ['file://' + f, os.path.join(self.s3_dir, tissue, os.path.basename(self.combat_tsv))])

    def teardown(self):
        pass

    def run_combat_script(self):
        p = 'T' if self.parametric else 'F'
        name = os.path.basename(self.combat_tsv)

        return textwrap.dedent("""
            source('combat.r')
            ComBat('{name}', 'batch.txt', skip=1, par.prior='{p}')
            """.format(**locals()))
