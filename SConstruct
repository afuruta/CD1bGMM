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

# download Adaptive and NN-distantce data
download = env.Command('data.tar.gz',
                       None,
                       'wget https://www.dropbox.com/s/3dl1acdhfqhrlih/data.tar.gz')
# unzip and delete archive
data = env.Command('data/',
                   download,
                   'tar -zxvf $SOURCE && rm $SOURCE')


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
                                                          'decision.pdf',
                                                          'roc.pdf']],
                             data,
                             'python bin/classify.py ' + outbase + ' $SOURCE/Adaptive/clinical_groups.tsv $SOURCE/NNdist/tmp.*wsplit*NNdistance.tsv')
