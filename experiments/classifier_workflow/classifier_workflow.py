import matplotlib

matplotlib.use('Agg')  # Implement linux safe backend for plotting

import argparse
import logging
import os
import pickle
import sys
from collections import defaultdict
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
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.linear_model import RandomizedLogisticRegression
from sklearn.model_selection import KFold, train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import scale
from sklearn.svm import SVC, LinearSVC
from sklearn.tree import DecisionTreeClassifier

from utils import cls, mkdir_p

log = logging.getLogger(__name__)


class ClassifierWorkflow():
    def __init__(self, X, y, folds, verbose):
        self.X = pd.read_csv(X, sep='\t')
        self.y = np.array([x.strip() for x in open(y).readlines()])
        self.folds = folds
        self.verbose = verbose
        self.features = self.X.columns
        self.step = None
        self.n_features = None
        self.ranks = dict()
        self.scores = defaultdict(list)
        self.output = 'output'
        mkdir_p(self.output)

        # Classifiers
        self.classifiers = {'GradientBoostingClassifier': GradientBoostingClassifier(verbose=self.verbose),
                            'AdaBoost': AdaBoostClassifier(),
                            'LinearSVC': LinearSVC(verbose=self.verbose),
                            'MLPClassifier': MLPClassifier(alpha=1, verbose=self.verbose),
                            'DecisionTreeClassifier': DecisionTreeClassifier(),
                            'GaussianNB': GaussianNB(),
                            'QuadraticDiscriminantAnalysis': QuadraticDiscriminantAnalysis(),
                            'LinearDiscriminantAnalysis': LinearDiscriminantAnalysis()}

        self.multicore_classifiers = {
            'RandomForest': RandomForestClassifier(n_jobs=-1, verbose=self.verbose),
            'ExtraTrees': ExtraTreesClassifier(n_jobs=-1, verbose=self.verbose),
            'RandomizedLogisticRegression': RandomizedLogisticRegression(n_jobs=-1, verbose=self.verbose),
            'GaussianProcessClassifier': GaussianProcessClassifier(n_jobs=-1),
            'KNeighborsClassifier': KNeighborsClassifier(n_jobs=-1)}

    def get_train_test(self):
        if self.folds == 1:
            X_train, X_test, y_train, y_test = train_test_split(self.X, self.y)
            yield scale(X_train), y_train, scale(X_test), y_test
        else:
            kf = KFold(self.folds)
            for train, test in kf.split(self.X, self.y):
                yield scale(np.array(self.X)[train, :]), self.y[train], scale(np.array(self.X)[test, :]), self.y[test]

    def feature_selection_and_fitting(self, X, y, n=0):
        # RFECV isn't an estimator, just used to estimate the optimal number of features
        self.step = np.ceil(len(self.features) * 1.0 / 20)
        rfecv = RFECV(estimator=SVC(kernel="linear"), step=self.step, cv=StratifiedKFold(2),
                      scoring='accuracy', n_jobs=-1, verbose=self.verbose)

        # Methods to estimate feature importance
        multi_methods = ['RandomForest', 'ExtraTrees']
        single_methods = ['GradientBoostingClassifier', 'AdaBoost', 'DecisionTreeClassifier',
                          'LinearSVC', 'LinearDiscriminantAnalysis']

        multi_clfs = [self.multicore_classifiers[x] for x in multi_methods] + [rfecv]
        single_clfs = [self.classifiers[x] for x in single_methods]

        # Train Classifiers
        log.info('Training classifiers for feature selection')
        clfs = [clf.fit(X, y) for clf in multi_clfs]
        # Create multiprocessing pool so we can parallelize fitting
        p = Pool(cpu_count())
        single_clfs = p.map(unwrap_self_return_fit,
                            zip([self] * len(single_clfs), [(clf, X, y) for clf in single_clfs]))
        p.close()
        p.join()
        clfs.extend(single_clfs)

        # Unpack fit classifiers and store their respective coefficients
        for method, clf in zip(multi_methods + ['RFECV'] + single_methods, clfs):
            if method == 'LinearSVC' or method == 'LinearDiscriminantAnalysis':
                self.ranks[method] = rank_to_dict(np.abs(clf.coef_)[0], self.features)
            elif method == 'RFECV':
                log.info("Optimal number of features : %d" % rfecv.n_features_)
                self.n_features = int(rfecv.n_features_)
            else:
                self.ranks[method] = rank_to_dict(clf.feature_importances_, self.features)

        # Compute mean of all methods for each feature
        log.info('Computing average score for each feature')
        r = {}
        for feature in self.features:
            log.debug('Computing average value for feature: {}'.format(feature))
            r[feature] = round(np.median([self.ranks[method][feature] for method in self.ranks.keys()]), 2)

        # Output Dataframe of features by importance score
        log.info('Creating Dataframe of Feature Importance')
        expr = pd.DataFrame()
        expr['feature'] = self.features
        for method in sorted(self.ranks.keys()):
            expr[method] = [self.ranks[method][f] for f in self.features]
        expr['Mean'] = [r[f] for f in self.features]
        expr.sort_values('Mean', ascending=False)
        expr.to_csv(os.path.join(self.output, 'feature-scores-{}.tsv'.format(n)), sep='\t', index=False)

        # Plot number of features VS. cross-validation scores
        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (nb of correct classifications)")
        plt.plot([(x + 1) * self.step for x in range(len(rfecv.grid_scores_))], rfecv.grid_scores_)
        plt.title('Optimal Number of Features for Classification')
        plt.savefig(os.path.join(self.output, 'features-vs-cv-scores-{}.png'.format(n)))

    def _return_fit(self, (clf, X, y)):
        clf.fit(X, y)
        return clf

    def score(self, X_train, y_train, X_test, y_test, n=0):
        best_clf = None
        best_name = None
        best_score = 0

        # Using information from the feature selection step we'll select top K features
        # And then retrain the classifiers on the reduced dataset
        log.info('Selecting top {} features based on RFECV'.format(self.n_features))
        skb = SelectKBest(k=self.n_features)
        skb.fit(X_train, y_train)
        mask = skb.get_support(indices=True)
        X_train = X_train[:, mask]
        X_test = X_test[:, mask]

        log.info('Retraining and Scoring Classifiers')
        for name in sorted(self.classifiers.keys()):
            # Retrain and score
            clf = self.classifiers[name]
            clf.fit(X_train, y_train)
            score = clf.score(X_test, y_test)
            self.scores[name].append(score)

            # Retain top-scoring CLF
            if score > best_score:
                best_score = score
                best_clf = clf
                best_name = name

        log.info('Seralizing top-scoring classifier: {}\tScore: {}'.format(best_name, best_score))
        with open(os.path.join(self.output, '{}-{}.pickle'.format(best_name, n)), 'w') as f:
            pickle.dump(best_clf, f)

    def save_scores(self):
        log.info('Saving scores for all classifiers')
        cols = ['Score-R{}'.format(x + 1) for x in xrange(self.folds)]
        with open(os.path.join(self.output, 'scores.tsv'), 'w') as f:
            f.write('Method\t{}\tAverage\n'.format('\t'.join(cols)))
            for k, v in sorted(self.scores.iteritems()):
                f.write('{}\t{}\t{}\n'.format(k, '\t'.join([str(x) for x in v]), np.mean(v)))


def unwrap_self_return_fit(arg, **kwarg):
    return ClassifierWorkflow._return_fit(*arg, **kwarg)


def rank_to_dict(ranks, names):
    minmax = MinMaxScaler()
    ranks = minmax.fit_transform(ranks.reshape(-1, 1))
    ranks = map(lambda x: round(x, 2), ranks)
    return dict(zip(names, ranks))


def main():
    """
    Classifier Workflow
     
    -- Steps --
    1. Dataframe and label vector are split into train / test set KFold n times
        2. Each train / test group is scaled to mean 0 unit variance 1
        3. Optimum # of features is determined and feature importance is scored by several classifiers
        4. Features are reduced using SelectKBest and Chi2 with the value from step 3.
        5. Classifiers are trained and scored. The best scoring classifier is serialized.
    6. Save scores from n runs to `scores.tsv`
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=False)

    parser.add_argument("-X", help="Path to dataframe TSV of samples by features. Ensure all categorical features "
                                   "have been One-Hot-Encoded. Will be read in by Pandas and expects column labels.")
    parser.add_argument("-y", type=str, help='Path to label vector for dataframe (X). One label per '
                                             'line. Will be read in directly with NumPy.')
    parser.add_argument('-n', type=int, default=3, help='Number of folds to perform workflow')
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        cls()
        parser.print_help()
        sys.exit(1)
    cls()

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    verbose = 1 if args.verbose else 0

    c = ClassifierWorkflow(args.X, args.y, args.n, verbose)
    for n, (X_train, y_train, X_test, y_test) in enumerate(c.get_train_test()):
        n = n + 1
        log.info('\nBeginning run for fold {} of {}'.format(n, args.n))
        c.feature_selection_and_fitting(X_train, y_train, n=n)
        c.score(X_train, y_train, X_test, y_test, n=n)
    c.save_scores()


if __name__ == '__main__':
    main()
