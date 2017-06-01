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
    return ClassifierWorkflow._return_fit(*arg, **kwarg)


class ClassifierWorkflow():
    def __init__(self, X, y, folds):
        self.X = pd.read_csv(X, sep='\t', index_col=0)
        self.y = np.array([x.strip() for x in open(y).readlines()])
        self.folds = folds
        self.features = X.columns
        self.samples = X.index
        self.n_features = None
        self.ranks = dict()
        self.scores = defaultdict()

        # Classifiers
        self.classifiers = {'GradientBoostingClassifier': GradientBoostingClassifier(verbose=1),
                            'AdaBoost': AdaBoostClassifier(),
                            'LinearSVC': LinearSVC(verbose=1),
                            'MLPClassifier': MLPClassifier(alpha=1, verbose=1),
                            'DecisionTreeClassifier': DecisionTreeClassifier(),
                            'GaussianNB': GaussianNB(),
                            'QuadraticDiscriminantAnalysis': QuadraticDiscriminantAnalysis(),
                            'LinearDiscriminantAnalysis': LinearDiscriminantAnalysis()}

        self.multicore_classifiers = {
            'RandomForest': RandomForestClassifier(n_jobs=-1, verbose=1),
            'ExtraTrees': ExtraTreesClassifier(n_jobs=-1, verbose=1),
            'RandomizedLogisticRegression': RandomizedLogisticRegression(n_jobs=-1, verbose=1),
            'GaussianProcessClassifier': GaussianProcessClassifier(n_jobs=-1),
            'KNeighborsClassifier': KNeighborsClassifier(n_jobs=-1)}

    def get_train_test(self):
        kf = KFold(self.folds)
        for train, test in kf.split(self.X, self.y):
            yield scale(self.X[train]), self.y[train], scale(self.X[test]), self.y[test]

    def feature_selection_and_fitting(self, X, y, n=0):
        # RFECV is a feature selection method not an estimator, so it's defined here
        rfecv = RFECV(estimator=SVC(kernel="linear"), step=1000, cv=StratifiedKFold(2),
                      scoring='accuracy', n_jobs=-1, verbose=1)

        # Methods to estimate feature importance
        multi_methods = ['RandomForest', 'ExtraTrees']
        single_methods = ['GradientBoostingClassifier', 'AdaBoost', 'DecisionTreeClassifier',
                          'LinearSVC', 'LinearDiscriminantAnalysis']

        multi_clfs = [self.multicore_classifiers[x] for x in multi_methods] + [rfecv]
        single_clfs = [self.classifiers[x] for x in single_methods]

        # Train Classifiers
        log.info('Training Classifiers')
        # Create multiprocessing pool so we can parallelize fitting
        p = Pool(cpu_count())
        clfs = [clf.fit(X, y) for clf in multi_clfs] + \
                p.map(unwrap_self_return_fit, zip([self] * len(single_clfs), [(clf, X, y) for clf in single_clfs]))

        # Unpack fit classifiers and store their respective coefficients
        for method, clf in zip(multi_methods + ['RFECV'] + single_methods, clfs):
            if method == 'LinearSVC' or method == 'LinearDiscriminantAnalysis':
                self.ranks[method] = rank_to_dict(np.abs(clf.coef_), self.features)
            elif method == 'RFE':
                log.info("Optimal number of features : %d" % rfecv.n_features_)
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
        expr.to_csv('feature-scores-{}.tsv'.format(n), sep='\t')

        # Plot number of features VS. cross-validation scores
        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (nb of correct classifications)")
        plt.plot([x * 1000 for x in range(len(rfecv.grid_scores_))], rfecv.grid_scores_)
        plt.title('Optimal Number of Features for Classification')
        plt.savefig('features-vs-cv-scores-{}.png'.format(n))

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
        skb = SelectKBest(chi2, k=self.n_features)
        skb.fit(X_train, y_train)
        mask = skb.get_support(indices=True)

        # Apply Mask
        X_train, y_train = X_train[mask], y_train[mask]
        X_test, y_test = X_test[mask], y_test[mask]

        log.info('\n\nRetraining and Scoring Classifiers')
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
                log.info('\tCurrent best classifier: {}'.format(name))

            log.info('Classifier: {}\tScore: {}'.format(name, score))

        log.info('Seralizing top-scoring classifier')
        with open('{}-{}.pickle'.format(best_name, n), 'w') as f:
            pickle.dump(best_clf, f)

    def save_scores(self):
        log.info('Saving scores for all classifiers')
        with open('scores.tsv', 'w') as f:
            f.write('Method\tScore-R1\tScore-R2\tScore-R3\tAverage\n')
            for k, v in sorted(self.scores.iteritems()):
                f.write('{}\t{}\t{}\n'.format(k, '\t'.join([str(x) for x in v]), np.mean(v)))


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
