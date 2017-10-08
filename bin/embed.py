#! /usr/bin/env python
# -*- coding: utf-8 -*-

from utils import *
import argparse, glob
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns; sns.set(style="white", color_codes=True)
sns.set_style('ticks')
import pandas as pd

def load_tsvs(tsvs, Vmap=None, downsample=None):
    '''load adaptive tsv files'''

    Mtb = [tsv_tags(tsv)['Bacterial Infections and Mycoses'] for tsv in tsvs]
    Mtb = scipy.array([x[0].lower() if len(x) == 1 else 'healthy' for x in Mtb])

    Vmap = getVmap(Vmap) if Vmap is not None else None
    v = DictVectorizer(sort=False)
    X = v.fit_transform(Vcdr3J(tsvs, Vmap, downsample=downsample))
    TCRs = scipy.array(v.feature_names_)

    return X, TCRs, Mtb

def main():

    scipy.seterr(divide='raise', over='raise',invalid='raise')

    parser = argparse.ArgumentParser('manifold embedding for TCRs')
    parser.add_argument('--out', type=str, default='embedding.pdf', help='output pdf file name')
    parser.add_argument('--jobs', type=int, default=1, help='number of concurent processes (careful of memory usage)')
    parser.add_argument('--perplexity', type=float, default=30., help='t-SNE perplexity')
    parser.add_argument('--Vmap', type=str, default=None, help='V map file for Adaptive to TCRdist nomenclature')
    parser.add_argument('tetramer', nargs='+', help='immunoseq tetramer sorted samples')
    args = parser.parse_args()

    print('generating sample/TCR matrix from immunoSEQ files:')
    X, TCRs, _ = load_tsvs(args.tetramer, Vmap=args.Vmap)
    print('done, {:,} samples X {:,} TCRs'.format(*X.shape), flush=True)

    y = []
    subject = []
    for tsv in args.tetramer:
        if '_re-sorted_to_99percent_Tetramer+' in tsv:
            new_y = 'resorted tet+'
        elif 'Tetramer+' in tsv:
            new_y = 'tet+'
        elif 'Tetramer-' in tsv:
            new_y = 'tet-'
        else:
            raise ValueError('unrecognized file type: ' + tsv)
        new_subject = os.path.basename(tsv).split('_')[0]
        print(new_subject, new_y)
        y.append(new_y)
        subject.append(new_subject)

    # filter out extremeley short CDR3s (they make TCRdist sad)
    long_enough = scipy.array([len(tcr[1]) > 5 for tcr in TCRs])
    X = X[:, long_enough]
    TCRs = TCRs[long_enough]

    where_are_NaNs = scipy.isnan(X.data)
    assert where_are_NaNs.sum() == 0

    distance_calculator = DistanceCalculator('human')

    distances = pairwise_distances(scipy.array([range(TCRs.shape[0])]).T, metric=index_dist, n_jobs=args.jobs, TCRs=TCRs, distance_calculator=distance_calculator)
    print('manifold learning... ', end='', flush=True)
    assert distances.min() >= 0

    X_manifold = TSNE(perplexity=args.perplexity, metric='precomputed', verbose=1).fit_transform(distances)
    df = pd.DataFrame()
    core_samples = []
    for i in range(X.shape[0]):
        df_dict = {}
        df_dict['TCR_indices'] = scipy.ravel(X[i, :].toarray()).nonzero()[0]
        df_dict['TSNE 1'] = X_manifold[df_dict['TCR_indices'], 0]
        df_dict['TSNE 2'] = X_manifold[df_dict['TCR_indices'], 1]
        df_dict['subject'] = subject[i]
        df_dict['sort'] = y[i]
        df_dict['abundance'] = scipy.ravel(X[i, df_dict['TCR_indices']].toarray())
        # make it fractional
        df_dict['abundance'] = df_dict['abundance']/df_dict['abundance'].sum()
        df = df.append(pd.DataFrame(df_dict))

    # add indicators for canonical motifs (Van Rhijn)
    canonicals = []
    for i in df['TCR_indices']:
        if (sum('TRAV1-2' in x for x in TCRs[i][0]) > 0 and len(TCRs[i][1]) == 13 and TCRs[i][2] == 'TCRAJ09'):
            canonicals.append('V1-2/J9')
        elif (sum('TRBV6-2' in x for x in TCRs[i][0]) > 0 and len(TCRs[i][1]) == 14 and TCRs[i][2] == 'TCRBJ02'):
            canonicals.append('V6-2/J2')
        else:
            canonicals.append(None)
    df['canonical'] =  canonicals

    plt.figure(figsize=(4, 4))
    for subject, color in zip(sorted(set(df['subject'])), ('blue', 'magenta', 'green', 'cyan')):
        sns.regplot(x='TSNE 1', y='TSNE 2', data=df.loc[(df['subject']==subject) & (df['sort'] != 'tet-')], fit_reg=False, scatter_kws={'alpha':.3}, color=color)
        sns.regplot(x='TSNE 1', y='TSNE 2', data=df.loc[(df['subject']==subject) & (df['sort'] != 'tet-') & df.canonical], fit_reg=False, scatter_kws={'alpha':.6, 'edgecolor':'k', 'linewidth':1, 'facecolors':'none'})
    plt.tight_layout()
    plt.savefig(args.out)


if __name__ == '__main__':
    main()
