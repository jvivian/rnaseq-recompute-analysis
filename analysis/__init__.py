import pandas as pd

tissues = ['adrenal',
           'bladder',
           'breast',
           'cervix',
           'colon',
           'esophagus',
           'kidney',
           'liver',
           'lung_adenocarcinoma',
           'lung_mesothelioma',
           'lung_squamous',
           'ovary',
           'pancreas',
           'prostate',
           'skin',
           'stomach',
           'testis',
           'thyroid',
           'uterus_carcinosarcoma',
           'uterus_endometrioid_carcinoma']

samples_df = pd.DataFrame()
samples_df['gtex'] = [125, 9, 178, 10, 140, 379, 27, 110, 287, 287, 287, 88, 165, 100, 556, 173, 165, 278, 78, 78]
samples_df['tcga_t'] = [77, 407, 1092, 303, 287, 181, 530, 369, 513, 87, 498, 420, 178, 494, 102, 413, 148, 504, 57, 180]
samples_df['tcga_n'] = [0, 19, 113, 3, 41, 13, 72, 50, 59, 0, 50, 0, 4, 51, 1, 36, 0, 59, 0, 23]
samples_df['tcga_2'] = [0, 0, 0, 0, 1, 0, 0, 2, 2, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 1]
samples_df['tcga_5'] = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0]
samples_df['tcga_6'] = [0, 0, 7, 2, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 366, 0, 0, 8, 0, 0]
samples_df.index = tissues
