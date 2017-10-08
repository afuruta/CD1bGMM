#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import os
from sklearn.metrics import log_loss, accuracy_score, f1_score, roc_curve, roc_auc_score
import pickle
from tcrdist.DistanceCalculator import DistanceTCR, DistanceCalculator
import scipy
from collections import defaultdict
from sklearn.feature_extraction import DictVectorizer

def index_dist(i, j, TCRs=None, distance_calculator=None):
        '''hack to work with sklearn and string data'''
        return TCRdist(TCRs[int(i[0])], TCRs[int(j[0])], distance_calculator=distance_calculator)

def TCRdist(key1, key2, distance_calculator=None):
    if distance_calculator is None:
        distance_calculator = DistanceCalculator('human')
    v_1, cdr3_1, _ = key1
    v_2, cdr3_2, _ = key2
    mismatch_error = ValueError('chains mismatch: {}, {}'.format(v_1, v_2))
    if 'TRB' in list(v_1)[0]:
        if 'TRB' not in list(v_2)[0]:
            raise mismatch_error
    elif 'TRA' in list(v_1)[0]:
        if 'TRA' not in list(v_2)[0]:
            raise mismatch_error
    x1 = DistanceTCR('human', [], v_1, '', cdr3_1)
    x2 = DistanceTCR('human', [], v_2, '', cdr3_2)
    result = distance_calculator.distance(x1, x2, chains='B')
    # if result == 0:
    #     assert tuple(key1) == tuple(key2)
    return result

def getVmap(mapfile):
    Vmap = {}
    with open(mapfile, 'r') as f:
        for line in f:
            phil_name, adaptive_name = line.rstrip().split()
            family, gene, allele = adaptive_name.split(',')
            if family not in Vmap:
                Vmap[family] = {gene:{allele:phil_name}}
            else:
                if gene not in Vmap[family]:
                    Vmap[family][gene] = {}
                Vmap[family][gene][allele] = phil_name
    return Vmap

def Vtranslate(family, genes, alleles, Vmap):
    if len(genes) > 1:
        assert alleles == None
        return frozenset(allele for gene in genes for allele in Vmap[family][gene].values())
    else:
        return frozenset(Vmap[family][genes[0]][allele] for allele in alleles)

def Vcdr3J(tsvs, Vmap=None, downsample=None):
    '''Adaptive tsvs to python dictionaries'''
    # theres a lot of repeated code here for purposes of optimization,
    # like not checking conditions on each iteration through the file lines
    # if Vmap is not None, get in the form Phil's code expects
    ntsvs = len(tsvs)
    for tsv_ct, tsv in enumerate(tsvs, 1):
        dat = {}
        uniques = 0
        template_total = 0
        badsample = False
        with open(tsv, 'r') as f:
            fields = f.readline().rstrip().split('\t')
            v_family_index = fields.index('v_family')
            v_gene_index = fields.index('v_gene')
            v_allele_index = fields.index('v_allele')
            j_family_index = fields.index('j_family')
            j_gene_index = fields.index('j_gene')
            j_allele_index = fields.index('j_allele')
            frame_type_index = fields.index('frame_type')
            amino_acid_index = fields.index('amino_acid')
            templates_index = fields.index('templates')
            if Vmap is not None:
                v_gene_ties_index = fields.index('v_gene_ties')
                v_allele_ties_index = fields.index('v_allele_ties')
            for line in f:
                fields = line.rstrip().split('\t')
                if fields[frame_type_index] != 'In': continue
                if not fields[templates_index].isdigit():
                    badsample = True
                    templates = scipy.nan
                else:
                    templates = int(fields[templates_index])
                    if templates == 0:
                        badsample = True
                        templates = scipy.nan
                if Vmap is not None:
                    family = fields[v_family_index]
                    # unresolved or pseudogenes
                    if family in ('', 'unresolved') or len(fields[amino_acid_index]) < 2: continue
                    genes = (fields[v_gene_index],) if fields[v_gene_index] not in ('', 'unresolved') else tuple(fields[v_gene_ties_index].split(','))
                    if len(genes) == 1:
                        alleles = (fields[v_allele_index],) if fields[v_allele_index] not in ('', 'unresolved') else tuple(fields[v_allele_ties_index].split(','))
                    else:
                        alleles = None
                    v_alleles = Vtranslate(family, genes, alleles, Vmap)
                    if 'TRBV22/OR9-2*01' in v_alleles or 'TRBV7-5*01' in v_alleles \
                        or 'TRBV7-5*02' in v_alleles or 'TRBV25/OR9-2*02' in v_alleles \
                        or 'TRBVA/OR9-2*01' in v_alleles or 'TRBV8-2*01' in v_alleles \
                        or 'TRBV22-1*01' in v_alleles or 'TRBV5-2*01' in v_alleles \
                        or 'TRDV1*01' in v_alleles or 'TRDV2*03' in v_alleles \
                        or 'TRAV8-5*01' in v_alleles or 'TRAV3*02' in v_alleles \
                        or 'TRAV33*01' in v_alleles or 'TRAV29/DV5*03' in v_alleles \
                        or 'TRDV3*01' in v_alleles or 'TRAV31*01' in v_alleles \
                        or 'TRAV37*01' in v_alleles or 'TRAV28*01' in v_alleles \
                        or 'TRAV32*01' in v_alleles or 'TRAV15*01' in v_alleles: continue
                    key = (v_alleles, fields[amino_acid_index], fields[j_family_index])
                else:
                    key = '|'.join([fields[v_family_index], fields[v_gene_index], fields[v_allele_index],
                                    fields[j_family_index], fields[j_gene_index], fields[j_allele_index],
                                    fields[amino_acid_index]])
                if key in dat:
                    dat[key] += templates
                else:
                    dat[key] = templates
                    uniques += 1
                template_total += templates

        print('    file {} of {}, '.format(tsv_ct, ntsvs) + os.path.basename(tsv) + ', uniques: {}, '.format(uniques) + ('WARNING: invalid template values' if badsample else 'templates: {}'.format(template_total)))

        if downsample is not None:
            yield dict(sorted(dat.items(), key=lambda x: x[1], reverse=True)[:downsample])
        else:
            yield dat


def tsv_tags(tsv):
    '''read Adaptive tsv sample tags'''
    with open(tsv, 'r') as f:
        fields = f.readline().rstrip().split('\t')
        d = defaultdict(list)
        for x in f.readline().rstrip().split('\t')[fields.index('sample_tags')].split(','):
            key, value = x.split(':')
            d[key].append(value)
        return d
