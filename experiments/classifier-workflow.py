import matplotlib
from collections import defaultdict

matplotlib.use('Agg')  # Implement linux safe backend for plotting

import argparse
import logging
import pickle
from multiprocessing import Pool
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.linear_model import RandomizedLogisticRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import scale
from sklearn.svm import SVC, LinearSVC
from sklearn.tree import DecisionTreeClassifier

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def unwrap_self_return_fit(arg, **kwarg):
    return Classifier._return_fit(*arg, **kwarg)


class Classifier():
    def __init__(self, X, y):
        self.X = pd.read_csv(X, sep='\t', index_col=0)
        self.y = np.array([x.strip() for x in open(y).readlines()])
        self.features = X.columns
        self.samples = X.index
        self.ranks = dict()

    def get_train_test(self, n=3):
        kf = KFold(n)
        for train, test in kf.split(self.X, self.y):
            yield scale(self.X[train]), self.y[train], scale(self.X[test]), self.y[test]

    def _return_fit(self, (clf, X, y, features)):
        clf.fit(X, y)
        return clf

    def feature_selection(self, X, y):

        # Create multiprocessing pool so we can parallelize fitting
        p = Pool(cpu_count())

        # Methods to use for feature selection
        multi_methods = [
            RandomForestClassifier(n_jobs=-1, verbose=1),
            ExtraTreesClassifier(n_jobs=-1, verbose=1),
            RFECV(estimator=SVC(kernel="linear"), step=1000,
                  cv=StratifiedKFold(2), scoring='accuracy', n_jobs=-1, verbose=1),
            RandomizedLogisticRegression(n_jobs=-1, verbose=1)
        ]

        # These methods can't use multicore so we'll parallelize them
        single_core_methods = [
            GradientBoostingClassifier(verbose=1),
            AdaBoostClassifier()
        ]

        log.info('Training Classifiers')
        clfs = [clf.fit(X, y) for clf in multi_methods] + \
               p.map(unwrap_self_return_fit, zip([self] * len(single_core_methods), single_core_methods))


        # Unpack fit classifiers and store their respective coefficients
        rf, et, rfecv, rlr, gbc, ada = clfs
        self.ranks["RF"] = rank_to_dict(np.abs(rf.feature_importances_), self.features)
        self.ranks["ETrees"] = rank_to_dict(np.abs(et.feature_importances_), self.features)
        self.ranks["RFECV"] = rank_to_dict(np.abs(rfecv.feature_importances_), self.features)
        self.ranks['RLR'] = rank_to_dict(np.abs(rlr.feature_importances_), self.features)
        self.ranks["GBC"] = rank_to_dict(np.abs(gbc.feature_importances_), self.features)
        self.ranks["ADA"] = rank_to_dict(np.abs(gbc.feature_importances_), self.features)
        log.info("Optimal number of features : %d" % rfecv.n_features_)

        log.debug('Computing average value for each feature')
        r = {}
        for name in self.features:
            r[name] = round(np.median([self.ranks[method][name] for method in self.ranks.keys()]), 2)

        methods = sorted(self.ranks.keys())
        self.ranks["Mean"] = r
        methods.append("Mean")

        # Output Dataframe of features by importance score
        log.info('Creating')
        expr = pd.DataFrame()
        expr['feature'] = self.features
        for method in methods:
            expr[method] = [self.ranks[method][f] for f in self.features]
        expr.to_csv('feature-scores.tsv', sep='\t')

        # Plot number of features VS. cross-validation scores
        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (nb of correct classifications)")
        plt.plot([x * 1000 for x in range(len(rfecv.grid_scores_))], rfecv.grid_scores_)
        plt.title('Optimal Number of Features for Classification')
        plt.savefig('features-vs-cv-scores.png')


def rank_to_dict(ranks, names, order=1):
    minmax = MinMaxScaler()
    ranks = minmax.fit_transform(order * np.array([ranks]).T).T[0]
    ranks = map(lambda x: round(x, 2), ranks)
    return dict(zip(names, ranks))


def main():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)

    parser.add_argument("-X", help="Path to dataframe TSV of samples by features. Ensure all categorical features "
                                   "have been One-Hot-Encoded or something similar. Will be read in by Pandas.")
    parser.add_argument("-y", require=True, type=str, help='Path to label vector for dataframe (X). One label per'
                                                           'line. Will be read in directly with NumPy.')

    args = parser.parse_args()
    c = Classifier(args.X, args.y)
    for X_train, y_train, X_test, y_test in c.get_train_test():
        c.features(X_train, y_train)
        break

if __name__ == '__main__':
    main()
