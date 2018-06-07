#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
generate output plots for CD1b-GMM paper
'''
from __future__ import print_function
import os
import sys
# import sconsutils
from warnings import warn
from SCons.Script import Environment, AddOption

# sconsutils

# Set up SCons environment
environ = os.environ.copy()
env = Environment(ENV=environ)

# Add bin dir to PATH
env.PrependENVPath('PATH', 'bin')

# download TCRdist resources
tcrdist_setup = env.Command(['bin/tcrdist/datasets',
                             'bin/tcrdist/db',
                             'bin/tcrdist/external',
                             'bin/tcrdist/testing_ref'],
                            None,
                            'cd bin/tcrdist && python setup.py')

# download/unzip Adaptive tsv and TCRdist NN-distantce data
immuneaccess_link = 'https://s3-us-west-2.amazonaws.com/publishedproject-supplements/seshadri-2018-journalofimmunology/seshadri-2018-journalofimmunology-data.tar.gz'
data = env.Command('data/',
                   None,
                   'mkdir data && wget {} -O data.tar.gz && tar -zxvf data.tar.gz && rm data.tar.gz'.format(immuneaccess_link))


# manifold embedding for CD1b-GMM tetramer sorted sequences
embed_TRA = env.Command('output/embed_TRA.pdf',
                        data,
                        'python bin/embed.py --out $TARGET --Vmap $SOURCE/Adaptive/Vmap.tsv $SOURCE/Adaptive/tetramer_sorts/CS*Tetramer+_TCRAD.tsv')
embed_TRB = env.Command('output/embed_TRB.pdf',
                        data,
                        'python bin/embed.py --out $TARGET --Vmap $SOURCE/Adaptive/Vmap.tsv $SOURCE/Adaptive/tetramer_sorts/CS*Tetramer+.tsv')
env.Depends(embed_TRA, tcrdist_setup)
env.Depends(embed_TRB, tcrdist_setup)

# classification results using GEM motifs and NN-distance
outbase = 'output/classify'
classification = env.Command([outbase + '.' + f for f in ['canonical.pdf',
                                                          'cdr3_length.pdf',
                                                          'meanNNdist.pdf',
                                                          'GEM-TCR_motifs.decision.pdf',
                                                          'TCRdist.decision.pdf',
                                                          'roc.pdf']],
                             data,
                             'python bin/classify.py ' + outbase + ' $SOURCE/Adaptive/clinical_groups.tsv $SOURCE/NNdist/tmp.*wsplit*NNdistance.tsv')
