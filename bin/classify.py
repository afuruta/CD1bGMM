#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import argparse, scipy
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.stats import ttest_ind
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import seaborn as sns; sns.set(style="white", color_codes=True)
sns.set_style('ticks')
import pandas as pd
from matplotlib.ticker import MaxNLocator

def parse_nndist(nndist_files, clinical_groups):
    '''parse raw NNdist files'''
    df_clinical_groups = pd.read_csv(clinical_groups, sep='\t', index_col=0)
    df_clinical_groups.loc[df_clinical_groups['class'] == 'negative'] = 'IGRA\nnegative'
    df_clinical_groups.loc[df_clinical_groups['class'] == 'latent'] = 'IGRA\npositive'
    df_clinical_groups.loc[df_clinical_groups['class'] == 'active'] = 'active\nTB'
    df = pd.DataFrame()
    for file in nndist_files:
        print('processing ' + file)
        if 'resort' in file:
            sorted_population = 'in vitro'
        elif 'tetpos' in file:
            sorted_population = 'ex vivo'
        else:
            raise ValueError('could not parse sorted population from file name: ' + file)
        if 'alpha' in file:
            locus = 'alpha'
        elif 'beta' in file:
            locus = 'beta'
        else:
            raise ValueError('could not parse locus from file name ' + file)
        df_new = pd.read_csv(file,
                             delim_whitespace=True,
                             usecols=(1, 5, 11, 12, 16, 19),
                             header=None,
                             names=('subject', 'nearest neighbor distance', 'V', 'CDR3', 'abundance', 'J'))
        df_new['sorted population'] = sorted_population
        df_new['locus'] = locus
        df_new['subject'] = df_new.subject.str.split('_').str.get(3)
        df = df.append(df_new, ignore_index=True)
    df['clinical group'] = list(df_clinical_groups.loc[list(df.subject), 'class'])
    df['CDR3 length'] = [len(cdr3) for cdr3 in df['CDR3']]

    return df


def main():

    parser = argparse.ArgumentParser('clinical classification using GEM TCR burden and NN-distance')
    parser.add_argument('outbase', type=str, default='classify', help='base name for output files')
    parser.add_argument('clinical_groups', type=str, help='clinical group file')
    parser.add_argument('nndist_files', type=str, nargs='+', help='NN-dist files')
    args = parser.parse_args()

    # parse raw nndist files
    df = parse_nndist(args.nndist_files, args.clinical_groups)

    # first some grouping for the Van Rhijn et al. motifs
    df_burden = df.loc[df['sorted population'] == 'ex vivo', ['subject', 'clinical group', 'locus', 'V', 'J', 'CDR3 length', 'sorted population', 'abundance']]
    grouped = df_burden.groupby(['subject', 'clinical group'])
    def aggregator(group):
        match_alpha = ((group.locus == 'alpha') &
                       group.V.str.contains('TRAV1-2\*') &
                       group.J.str.contains('TRAJ9\*') &
                       (group['CDR3 length'] == 13))
        match_beta = ((group.locus == 'beta') &
                      group.V.str.contains('TRBV6-2\*') &
                      group.J.str.contains('TRBJ2-') &
                      (group['CDR3 length'] == 14))

        total_alpha = sum(group.locus == 'alpha')
        total_beta  = sum(group.locus == 'beta')
        return pd.Series({r'GEM-TCR$\alpha$' '\nmotif frequency':match_alpha.sum()/total_alpha,
                          r'GEM-TCR$\beta$' '\nmotif frequency':match_beta.sum()/total_beta})
    df_burden = grouped.apply(aggregator).reset_index()
    df_burden_rows_1052 = df_burden.subject.isin(['TB-1052', 'TB-1052-2M', 'TB-1052-6M'])

    print('alpha GEM-TCR motif significance:')
    print('    latent Vs negative: {}'.format(ttest_ind(df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\nnegative'), r'GEM-TCR$\alpha$' '\nmotif frequency'],
                                                    df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\npositive'), r'GEM-TCR$\alpha$' '\nmotif frequency'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\nnegative'), r'GEM-TCR$\alpha$' '\nmotif frequency'],
                                                    df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'active\nTB'), r'GEM-TCR$\alpha$' '\nmotif frequency'])[1]))
    print('beta GEM-TCR motif significance:')
    print('    latent Vs negative: {}'.format(ttest_ind(df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\nnegative'), r'GEM-TCR$\beta$' '\nmotif frequency'],
                                                    df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\npositive'), r'GEM-TCR$\beta$' '\nmotif frequency'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'IGRA\nnegative'), r'GEM-TCR$\beta$' '\nmotif frequency'],
                                                    df_burden.loc[~df_burden_rows_1052 & (df_burden['clinical group'] == 'active\nTB'), r'GEM-TCR$\beta$' '\nmotif frequency'])[1]))


    plt.figure(figsize=(4,4))
    plt.subplot(2, 1, 1)
    ax = sns.boxplot(x='clinical group', y=r'GEM-TCR$\alpha$' '\nmotif frequency', data=df_burden[~df_burden_rows_1052],
                     order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='clinical group', y=r'GEM-TCR$\alpha$' '\nmotif frequency', data=df_burden[~df_burden_rows_1052],
                  order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), dodge=True, edgecolor='k', linewidth=1)
    plt.legend(loc=(.01, .4), frameon=True)
    plt.ylim([0, None])
    sns.despine()
    plt.tight_layout()
    plt.subplot(2, 1, 2)
    ax = sns.boxplot(x='clinical group', y=r'GEM-TCR$\beta$' '\nmotif frequency', data=df_burden[~df_burden_rows_1052],
                     order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='clinical group', y=r'GEM-TCR$\beta$' '\nmotif frequency', data=df_burden[~df_burden_rows_1052],
                  order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), dodge=True, edgecolor='k', linewidth=1)
    plt.ylim([0, None])
    sns.despine()
    plt.tight_layout()
    plt.savefig(args.outbase + '.canonical.pdf')


    grouped = df.groupby(['subject', 'clinical group', 'locus', 'sorted population'])
    # aggregation function for mean, and total abundance
    def aggregator(group):
        total_abundance = group['abundance'].sum()
        mean_nndist = group['nearest neighbor distance'].mean()
        mean_cdr3length = group['CDR3 length'].mean()
        return pd.Series({'CD1b-GMM repertoire distance':mean_nndist,
                          'mean CDR3 length':mean_cdr3length,
                          'total abundance':total_abundance})
    df_mean = grouped.apply(aggregator).reset_index()

    # need to keep track of 1052 for exclusion
    df_rows_1052 = df.subject.isin(['TB-1052', 'TB-1052-2M', 'TB-1052-6M'])
    df_mean_rows_1052 = df_mean.subject.isin(['TB-1052', 'TB-1052-2M', 'TB-1052-6M'])

    # mean cdr3 length distributions
    plt.figure(figsize=(4,5))
    plt.subplot(2, 1, 1)
    ax = sns.boxplot(x='clinical group', y='mean CDR3 length', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population']=='ex vivo')],
                     order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='clinical group', hue='clinical group', y='mean CDR3 length', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population']=='ex vivo')],
                  order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), edgecolor='k', linewidth=1)
    ax.legend_.remove()
    plt.ylabel(r'TCR$\alpha$' '\nmean CDR3 length')
    #plt.ylim([0, None])
    sns.despine()
    plt.tight_layout()
    plt.subplot(2, 1, 2)
    ax = sns.boxplot(x='clinical group', y='mean CDR3 length', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population']=='ex vivo')],
                     order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='clinical group', hue='clinical group', y='mean CDR3 length', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population']=='ex vivo')],
                  order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), edgecolor='k', linewidth=1)
    ax.legend_.remove()
    plt.ylabel(r'TCR$\beta$' '\nmean CDR3 length')
    #plt.ylim([0, None])
    sns.despine()
    plt.tight_layout()
    plt.savefig(args.outbase + '.cdr3_length.pdf', bbox_inches='tight')

    # mean NNdists
    print('alpha TCRdist distance significance:')
    print('ex vivo')
    print('    latent Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\npositive'), 'CD1b-GMM repertoire distance'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'active\nTB')    , 'CD1b-GMM repertoire distance'])[1]))
    print('in vitro')
    print('    latent Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\npositive'), 'CD1b-GMM repertoire distance'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'ex vivo') & (df_mean['clinical group'] == 'active\nTB')    , 'CD1b-GMM repertoire distance'])[1]))
    print('beta TCRdist distance significance:')
    print('ex vivo')
    print('    latent Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\npositive'), 'CD1b-GMM repertoire distance'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='alpha') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'active\nTB')    , 'CD1b-GMM repertoire distance'])[1]))
    print('in vitro')
    print('    latent Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\npositive'), 'CD1b-GMM repertoire distance'])[1]))
    print('    active Vs negative: {}'.format(ttest_ind(df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'IGRA\nnegative'), 'CD1b-GMM repertoire distance'],
                                                        df_mean.loc[(~df_mean_rows_1052) & (df_mean.locus=='beta') & (df_mean['sorted population'] == 'in vitro') & (df_mean['clinical group'] == 'active\nTB')    , 'CD1b-GMM repertoire distance'])[1]))


    plt.figure(figsize=(5,3))
    plt.subplot(1, 2, 1)
    ax = sns.boxplot(x='sorted population', hue='clinical group', y='CD1b-GMM repertoire distance', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='alpha')],
                hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='sorted population', hue='clinical group', y='CD1b-GMM repertoire distance', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='alpha')],
                  hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), dodge=True, edgecolor='k', linewidth=1)
    ax.legend_.remove()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    sns.despine()
    plt.tight_layout()
    plt.subplot(1, 2, 2)
    ax = sns.boxplot(x='sorted population', hue='clinical group', y='CD1b-GMM repertoire distance', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='beta')],
                hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    sns.swarmplot(x='sorted population', hue='clinical group', y='CD1b-GMM repertoire distance', data=df_mean[(~df_mean_rows_1052) & (df_mean.locus=='beta')],
                hue_order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), dodge=True, edgecolor='k', linewidth=1)
    ax.legend_.remove()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    sns.despine()
    plt.tight_layout()
    plt.savefig(args.outbase + '.meanNNdist.pdf')

    df2 = pd.DataFrame(columns=('subject', 'clinical group', 'total alpha abundance', 'total beta abundance', r'mean TCR$\alpha$ CDR3 length', r'mean TCR$\beta$ CDR3 length', r'GEM-TCR$\alpha$' '\nmotif frequency', r'GEM-TCR$\beta$' '\nmotif frequency', 'x1', 'x2', 'x3', 'x4'))
    for i, (name, group) in enumerate(df_mean.groupby(('subject', 'clinical group'))):
        df2.loc[i] = list(name) + \
                     [group['total abundance'][group.locus == 'alpha'].iloc[0],
                      group['total abundance'][group.locus == 'beta'].iloc[0],
                      group['mean CDR3 length'][group.locus == 'alpha'].iloc[0],
                      group['mean CDR3 length'][group.locus == 'beta'].iloc[0],
                      df_burden.loc[df_burden['subject'] == name[0], r'GEM-TCR$\alpha$' '\nmotif frequency'].iloc[0],
                      df_burden.loc[df_burden['subject'] == name[0], r'GEM-TCR$\beta$' '\nmotif frequency'].iloc[0],
                      group['CD1b-GMM repertoire distance'][(group['sorted population'] == 'ex vivo')  & (group.locus == 'alpha')].iloc[0],
                      group['CD1b-GMM repertoire distance'][(group['sorted population'] == 'in vitro') & (group.locus == 'alpha')].iloc[0],
                      group['CD1b-GMM repertoire distance'][(group['sorted population'] == 'ex vivo')  & (group.locus == 'beta')].iloc[0],
                      group['CD1b-GMM repertoire distance'][(group['sorted population'] == 'in vitro') & (group.locus == 'beta')].iloc[0]]

    df2_rows_1052 = df2.subject.isin(['TB-1052', 'TB-1052-2M', 'TB-1052-6M'])
    X_canonical = df2.as_matrix(columns=(r'GEM-TCR$\alpha$' '\nmotif frequency', r'GEM-TCR$\beta$' '\nmotif frequency'))
    X_nndist = df2.as_matrix(columns=('x1', 'x2', 'x3', 'x4'))
    y = df2['clinical group'].as_matrix()
    # fit on only active vs negative, and no 1052

    df_roc = pd.DataFrame()
    for method, X in zip(['GEM-TCR motifs', 'TCRdist'], [X_canonical, X_nndist]):
        clf = LinearDiscriminantAnalysis()
        clf.fit(X[~df2_rows_1052, :], y[~df2_rows_1052])
        df2.loc[df2_rows_1052, 'active TB probability'] = clf.predict_proba(X[df2_rows_1052, :])[:, list(clf.classes_).index('active\nTB')]
        df2.loc[df2_rows_1052, 'prediction'] = clf.predict(X[df2_rows_1052, :])
        cv_probs = cross_val_predict(clf,
                                     X[~df2_rows_1052, :],
                                     y[~df2_rows_1052],
                                     method='predict_proba',
                                     cv=10)
        df2.loc[~df2_rows_1052, 'active TB probability'] = cv_probs[:, 2]
        df2.loc[~df2_rows_1052, 'prediction'] = cross_val_predict(clf,
                                                                  X[~df2_rows_1052, :],
                                                                  y[~df2_rows_1052],
                                                                  cv=10)
        print(method)
        print(confusion_matrix(y[~df2_rows_1052], df2.loc[~df2_rows_1052, 'prediction'].as_matrix(), labels=['IGRA\nnegative', 'IGRA\npositive', 'active\nTB']))


        df2.loc[df2.subject == 'TB-1052', 'time point'] = 'pre-tx'
        df2.loc[df2.subject == 'TB-1052-2M', 'time point'] = '2 mo'
        df2.loc[df2.subject == 'TB-1052-6M', 'time point'] = '6 mo'

        plt.figure(figsize=(3, 3))
        ax = sns.boxplot(data=df2[~df2_rows_1052], x='clinical group', y='active TB probability', order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), color='white')
        plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
        plt.setp(ax.lines, color='k')
        sns.swarmplot(data=df2[~df2_rows_1052], x='clinical group', y='active TB probability', order=('IGRA\nnegative', 'IGRA\npositive', 'active\nTB'), palette=('gray', 'orange', 'red'), edgecolor='k', linewidth=1, clip_on=False)
        plt.axhline(y=.5, ls='--', color='k', lw=1)
        plt.ylim([0, 1])
        sns.despine()
        plt.savefig('{}.{}.decision.pdf'.format(args.outbase, method).replace(' ', '_'), bbox_inches='tight')

        # cross validated roc (negative and active only)
        fpr, tpr, _ = roc_curve(y[~df2_rows_1052] == 'active\nTB',
                                df2.loc[~df2_rows_1052, 'active TB probability'])
        fpr = scipy.insert(fpr, 0, 0.)
        tpr = scipy.insert(tpr, 0, 0.)
        df_roc = df_roc.append(pd.DataFrame({'false positive rate':fpr, 'true positive rate':tpr, 'method':'{}\n(AUC = {:.2})'.format(method, auc(fpr, tpr))}))

    plt.figure()
    g = sns.FacetGrid(data=df_roc, hue='method', size=3)
    g = g.map(plt.plot, 'false positive rate', 'true positive rate', ls='-', lw=4, alpha=.8, clip_on=False)
    plt.legend(loc='lower right', fontsize='xx-small')
    g.set(xticks=(0., .2, .4, .6, .8, 1), yticks=(0., .2, .4, .6, .8, 1))
    plt.plot((0, 1.), (0, 1), ls='--', c='black', lw=.5, zorder=0)
    plt.xlim([0, 1.02])
    plt.ylim([0, 1.02])
    plt.tight_layout()
    plt.savefig(args.outbase + '.roc.pdf')


if __name__ == '__main__':
    main()
