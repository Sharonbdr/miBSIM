import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import itertools as it
import pandas as pd
from scipy import stats
import copy

# cnn utils
def rev_comp(seq):
    """ Get reverse complement """

    match_dict = {'A': 'T',
                  'T': 'A',
                  'C': 'G',
                  'G': 'C'}

    return ''.join([match_dict[x] for x in seq][::-1])


def get_centered_stype(site8, seq):
    if site8 == seq[2:-2]:
        return '8mer'
    elif site8[:-1] == seq[2:-3]:
        return '7mer-m8'
    elif site8[1:] == seq[3:-2]:
        return '7mer-a1'
    elif site8[1:-1] == seq[3:-3]:
        return '6mer'
    elif site8[:-2] == seq[2:-4]:
        return '6mer-m8'
    elif site8[2:] == seq[4:-2]:
        return '6mer-a1'
    else:
        return 'no site'


def get_best_stype(site8, seq):
    if site8 in seq:
        return '8mer'
    elif site8[:-1] in seq:
        return '7mer-m8'
    elif site8[1:] in seq:
        return '7mer-a1'
    elif site8[1:-1] in seq:
        return '6mer'
    elif site8[:-2] in seq:
        return '6mer-m8'
    elif site8[2:] in seq:
        return '6mer-a1'
    else:
        return 'no site'


def one_hot_encode(seq):
    if len(seq) == 0:
        return []
    """ 1-hot encode ATCG sequence """
    nt_dict = {
        'A': 0,
        'T': 1,
        'C': 2,
        'G': 3,
        'X': 4
    }
    targets = np.ones([5, 4]) / 4.0
    targets[:4, :] = np.eye(4)
    seq = [nt_dict[nt] for nt in seq]
    return list(targets[seq].flatten())


def generate_random_seq(length):
    """Generate random sequence"""
    nts = ['A', 'T', 'C', 'G']
    seq = np.random.choice(nts, size=length, replace=True)
    return ''.join(seq)


def get_target_no_match(mirna_sequence, length):
    """Given a miRNA sequence, return a random target sequence without 4 nt of contiguous pairing"""
    rc = rev_comp(mirna_sequence[1:8]) + 'A'
    off_limits = [rc[ix:ix + 4] for ix in range(5)]
    while True:
        target = generate_random_seq(length)
        keep = True
        for subseq in off_limits:
            if subseq in target:
                keep = False
                break

        if keep:
            return target

def get_mir_no_match(site_sequence, length):
    """Given a target sequence, return a random miRNA sequence without 4 nt of contiguous pairing"""
    rc = rev_comp(site_sequence[2:-3])
    off_limits = [rc[ix:ix + 4] for ix in range(4)]
    while True:
        target = generate_random_seq(length)
        keep = True
        for subseq in off_limits:
            if subseq in target:
                keep = False
                break

        if keep:
            return target


# rnaplfold utils
def generate_12mers(site8, only_canon):
    mers = []
    if only_canon:
        all_6mers = ["".join(kmer) for kmer in list(it.product(["A", "C", "G", "T"], repeat=6))]
        for i in range(3):
            subseq = site8[i:i + 6]
            mers += [x[:i + 2] + subseq + x[i + 2:] for x in all_6mers]

    else:
        all_8mers = ["".join(kmer) for kmer in list(it.product(["A", "C", "G", "T"], repeat=8))]
        for i in range(5):
            subseq = site8[i:i + 4]
            mers += [x[:i + 2] + subseq + x[i + 2:] for x in all_8mers]

    mers = list(set(mers))
    return sorted(mers)


# get_features utils
def site_to_ints(site):
    nt_dict = {
        'A': 0,
        'T': 1,
        'C': 2,
        'G': 3,
        'X': -1
    }

    return [nt_dict[x] for x in site]

def mir_site_pair_to_ints(mir, site):
    nt_dict = {
        'A': 0,
        'T': 1,
        'C': 2,
        'G': 3,
    }

    ints = []
    for nt1 in mir:
        if nt1 == 'X':
            ints += [-1] * len(site)
            continue
        for nt2 in site:
            if nt2 == 'X':
                ints.append(-1)

            else:
                ints.append((nt_dict[nt1] * 4) + nt_dict[nt2])

    return ints

# def get_nsites(features):
#     nsites = features.reset_index()
#     nsites['nsites'] = 1
#     nsites = nsites.groupby(['transcript','mir']).agg({'nsites': np.sum})
#     return nsites


def sigmoid(vals):
    return 1.0 / (1.0 + np.exp(-1 * vals))


# def norm_matrix(mat):
#     means = np.mean(mat, axis=1).reshape([-1, 1])
#     return mat - means


# def get_r2_unnormed(preds, labels):
#     preds_normed = norm_matrix(preds)
#     labels_normed = norm_matrix(labels)
#     return calc_r2(preds_normed, labels_normed)

# biochem model utils

def expand_features_4D(transcripts, mirs, max_nsites, feature_list, feature_df):
    features_4D = np.zeros([len(transcripts), len(mirs), max_nsites, len(feature_list)])
    for ix, transcript in enumerate(transcripts):
        for iy, mir in enumerate(mirs):
            try:
                temp = feature_df.loc[(transcript, mir)]
                nsites = len(temp)
                features_4D[ix, iy, :nsites, :] = temp[feature_list].values
            except KeyError:
                continue

    mask = ((np.abs(np.sum(features_4D, axis=3))) != 0).astype(int)

    return features_4D, mask


def split_vals(vals_4D, zero_indices):
    """
    Given a 4D matrix with ka values and features, split into ka_vals (3D), features (4D), and nosite_features (4D)
    """

    ka_vals_3D = vals_4D[:, :, :, 0]
    features_4D = vals_4D[:, :, :, 1:]
    nosite_features_4D = copy.copy(features_4D)
    for ix in zero_indices:
        nosite_features_4D[:, :, :, ix - 1] = 0

    return ka_vals_3D, features_4D, nosite_features_4D