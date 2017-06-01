import shutil

from experiments.classifier_workflow.classifier_workflow import ClassifierWorkflow

def test_workflow(tmpdir):
    x = 'test-X.tsv'
    y = 'test-y.tsv'
    n = 3
    verbose = 0
    c = ClassifierWorkflow(x, y, n, verbose)
    for n, (X_train, y_train, X_test, y_test) in enumerate(c.get_train_test()):
        n = n + 1
        c.feature_selection_and_fitting(X_train, y_train, n=n)
        c.score(X_train, y_train, X_test, y_test, n=n)
    c.save_scores()
    shutil.rmtree('output')