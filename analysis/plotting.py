import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_overlap(d1, d2=None, tissues=None, output_path=None):
    """
    Plot overlap of values between dictionaries, or against itself.

    :param dict d1: Tissues as keys
    :param dict d2: Optional - A second dataframe to compare against
    :param list tissues: Optional - A list of tissues to use instead of keys from d1
    :param str output_path: Optional - Path for output. Defaults to 'plot_overlap.png' in cwd
    """
    overlap = []
    tissues = tissues if tissues else d1.keys()
    d2 = d2 if d2 else d1
    n = 0

    for tissue in tissues:
        n = len(d1[tissue]) if len(d1[tissue]) > 0 else 1
        overlap.append([float(len(set(d1[tissue]).intersection(set(d2[x])))) / n for x in tissues])

    overlap = pd.DataFrame(overlap)
    overlap.index = label_fix(tissues)
    overlap.columns = label_fix(tissues)

    fig, ax = plt.subplots(figsize=[12, 10])
    sns.heatmap(overlap, ax=ax, cmap='Blues')
    plt.title('Top {} Overlap Between Tissues'.format(n))
    if output_path:
        plt.savefig(output_path, dpi=300, format='png')


def plot_barplot_tissues(tissue_dict, output_path=None, title='Number of Significant Genes'):
    """
    Barplot for values in a dictionary. Keys=tissues, values=int

    :param dict tissue_dict: Dictionary of tissues with some associated value.
    """
    f, ax = plt.subplots(figsize=(16, 8))
    x, y = zip_sort(tissue_dict.keys(), [len(tissue_dict[x]) for x in tissue_dict.keys()])
    sns.barplot(label_fix(x), y, ax=ax)
    plt.title(title)
    if output_path:
        plt.savefig(output_path)


def plot_barplot_intersection(td1, td2, output_path=None):
    """
    Barplot for intersection of two dictionary's values. Keys between dictionaries should be identical.
    :param dict td1:
    :param dict td2:
    :param output_path:
    :return:
    """
    inter = []
    for t in sorted(td1.keys()):
        n = len(td1[t]) if len(td1[t]) > 0 else 1
        inter.append(len(td1[t].intersection(td2[t]))*1.0 / n)

    f, ax = plt.subplots(figsize=(16,8))
    sns.barplot(label_fix(sorted(td1.keys())), inter, ax=ax)
    if output_path:
        plt.savefig(output_path)


def zip_sort(x, y):
    vals = sorted(zip(x, y), key=lambda z: z[1])
    return zip(*vals)


def label_fix(l):
    return [x.replace('_', '\n').capitalize() for x in l]
