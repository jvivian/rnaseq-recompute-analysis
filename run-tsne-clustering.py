# coding: utf-8
"""
John Vivian
November 2016
"""
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import argparse
import logging
import os
import pickle
import shutil
import textwrap

import matplotlib.pyplot as plt
from tqdm import tqdm

from tissue_clustering.tsne_clustering import create_classification_vector
from tissue_clustering.tsne_clustering import create_combined_df
from tissue_clustering.tsne_clustering import find_protein_coding_genes
from tissue_clustering.tsne_clustering import plot_dimensionality_reduction
from tissue_clustering.tsne_clustering import run_tsne
from tissue_clustering.tsne_clustering import split_tcga_tumor_normal
from utility_functions import mkdir_p

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

tissues = {
    'adrenal': ['Adrenal_Gland.tsv',
                'Adrenocortical_carcinoma.tsv',
                'Pheochromocytoma_and_Paraganglioma.tsv'],
    'bladder': ['Bladder.tsv',
                'Bladder_Urothelial_Carcinoma.tsv'],
    'brain': ['Brain_Amygdala.tsv',
              'Brain_Anterior_cingulate_cortex_(BA24).tsv',
              'Brain_Caudate_(basal_ganglia).tsv',
              'Brain_Cerebellar_Hemisphere.tsv',
              'Brain_Cerebellum.tsv',
              'Brain_Cortex.tsv',
              'Brain_Frontal_Cortex_(BA9).tsv',
              'Brain_Hippocampus.tsv',
              'Brain_Hypothalamus.tsv',
              'Brain_Lower_Grade_Glioma.tsv',
              'Brain_Nucleus_accumbens_(basal_ganglia).tsv',
              'Brain_Putamen_(basal_ganglia).tsv',
              'Brain_Spinal_cord_(cervical_c_1).tsv',
              'Brain_Substantia_nigra.tsv',
              'Glioblastoma_multiforme.tsv', ],
    'breast': ['Breast_Mammary_Tissue.tsv',
               'Breast_invasive_carcinoma.tsv'],
    'cervix': ['Cervical_squamous_cell_carcinoma_and_endocervical_adenocarcinoma.tsv',
               'Cervix_Endocervix.tsv',
               'Cervix_Ectocervix.tsv'],
    'colon': ['Colon_Sigmoid.tsv',
              'Colon_Transverse.tsv',
              'Colon_adenocarcinoma.tsv'],
    'esophagus': ['Esophageal_carcinoma.tsv',
                  'Esophagus_Gastroesophageal_Junction.tsv',
                  'Esophagus_Mucosa.tsv',
                  'Esophagus_Muscularis.tsv'],
    'kidney': ['Kidney_Chromophobe.tsv',
               'Kidney_Cortex.tsv',
               'Kidney_renal_clear_cell_carcinoma.tsv',
               'Kidney_renal_papillary_cell_carcinoma.tsv'],
    'liver': ['Liver.tsv',
              'Liver_hepatocellular_carcinoma.tsv'],
    'lung': ['Lung.tsv',
             'Lung_adenocarcinoma.tsv',
             'Lung_squamous_cell_carcinoma.tsv',
             'Mesothelioma.tsv'],
    'ovary': ['Ovarian_serous_cystadenocarcinoma.tsv',
              'Ovary.tsv'],
    'pancreas': ['Pancreas.tsv',
                 'Pancreatic_adenocarcinoma.tsv'],
    'prostate': ['Prostate.tsv',
                 'Prostate_adenocarcinoma.tsv'],
    'skin': ['Skin_Cutaneous_Melanoma.tsv',
             'Skin_Not_Sun_Exposed_(Suprapubic).tsv',
             'Skin_Sun_Exposed_(Lower_leg).tsv'],
    'stomach': ['Stomach.tsv',
                'Stomach_adenocarcinoma.tsv'],
    'testis': ['Testis.tsv',
               'Testicular_Germ_Cell_Tumors.tsv'],
    'thyroid': ['Thyroid.tsv',
                'Thyroid_carcinoma.tsv'],
    'uterus': ['Uterine_Carcinosarcoma.tsv',
               'Uterine_Corpus_Endometrioid_Carcinoma.tsv',
               'Uterus.tsv']}


def main():
    """
    First run create_project.py, located in the root diretctory of this project.

    Runs t-SNE clustering for a set of tissues
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna-seq-directory', type=str, required=True,
                        help='Full path to rna-seq_analysis directory, created by running create_project.py')
    params = parser.parse_args()

    log.info(textwrap.dedent("""
████████╗   ███████╗███╗   ██╗███████╗    ██████╗ ██╗      ██████╗ ████████╗███████╗
╚══██╔══╝   ██╔════╝████╗  ██║██╔════╝    ██╔══██╗██║     ██╔═══██╗╚══██╔══╝██╔════╝
   ██║█████╗███████╗██╔██╗ ██║█████╗      ██████╔╝██║     ██║   ██║   ██║   ███████╗
   ██║╚════╝╚════██║██║╚██╗██║██╔══╝      ██╔═══╝ ██║     ██║   ██║   ██║   ╚════██║
   ██║      ███████║██║ ╚████║███████╗    ██║     ███████╗╚██████╔╝   ██║   ███████║
   ╚═╝      ╚══════╝╚═╝  ╚═══╝╚══════╝    ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   ╚══════╝
    """))

    log.info('Creating experiment diretories')
    root_dir = params.rna_seq_directory
    experiment_dir = os.path.join(root_dir, 'experiments/tsne-clustering')
    plot_dir = os.path.join(experiment_dir, 'plots')
    [mkdir_p(os.path.join(experiment_dir, x)) for x in tissues.keys()]
    mkdir_p(plot_dir)

    log.info('Copying over requisite dataframes')
    tissue_dataframes = os.path.join(root_dir, 'data/tissue-dataframes')
    for tissue in tqdm(tissues):
        for tsv in tissues[tissue]:
            if not os.path.exists(os.path.join(experiment_dir, tissue, tsv)):
                shutil.copy(os.path.join(tissue_dataframes, tsv), os.path.join(experiment_dir, tissue, tsv))

    log.info('Collecting t-SNE data for each tissue')
    tsne_output = os.path.join(experiment_dir, 'tsne.pickle')
    if os.path.exists(tsne_output):
        log.info('Pickle file found, loading: ' + tsne_output)
        tsne = pickle.load(open(tsne_output, 'rb'))
    else:
        log.info('Curating list of protein coding genes')
        pc_genes = find_protein_coding_genes(os.path.join(root_dir, 'metadata/gencode.v23.annotation.gtf'))
        tsne = {}
        for tissue in tqdm(tissues):
            files = split_tcga_tumor_normal(os.path.join(experiment_dir, tissue))
            vector, label = create_classification_vector(files)
            df = create_combined_df(files, pc_genes)
            x = run_tsne(df)
            tsne[tissue] = (x, files, label)
        with open(tsne_output, 'w') as f:
            pickle.dump(tsne, f)

    log.info('Creating one large subplot')
    f, axes = plt.subplots(len(tissues), figsize=(8, 64))
    cm = plt.get_cmap('Accent')
    for i, tissue in enumerate(sorted(tsne)):
        matrix, files, label = tsne[tissue]
        color_set = (cm(1. * z / len(files)) for z in xrange(len(files)))
        names = [os.path.basename(x).split('.tsv')[0] for x in files]
        length = [x for x in xrange(len(files))]
        for color, l, target_name in zip(color_set, length, names):
            axes[i].scatter(matrix[label == l, 0], matrix[label == l, 1], alpha=0.5, color=color, label=target_name)
        axes[i].legend(loc='best', fontsize=6)
        axes[i].set_title(tissue)
    plt.savefig(os.path.join(plot_dir, 'tsne-plots.pdf'), format='pdf')

    log.info('Creating per-tissue plots')
    for tissue in tsne:
        f, ax = plt.subplots()
        x, files, label = tsne[tissue]
        plot_dimensionality_reduction(ax, x, files, label, title=tissue)
        f.savefig(os.path.join(plot_dir, tissue + '.pdf'), format='pdf')


if __name__ == '__main__':
    main()
