import os

import numpy as np
import pandas as pd

import utils

def iter_res(mirnames, site8s, nbins, infile_seqs, infile_bg, num_bg, outfile):
    for mirname, site8 in zip(mirnames, site8s):
        SA_bg = []
        for ix in range(nbins):
            seqs = pd.read_csv(os.path.join(infile_seqs.replace('MIR', mirname).replace('IX', str(ix))),
                               sep='\t')
            temp = pd.read_csv(os.path.join(infile_bg.replace('MIR', mirname).replace('IX', str(ix))), sep='\t',
                               header=None)
            temp.columns = ['12mer', 'p', 'logp']
            temp['count'] = 1
            temp = temp.groupby('12mer').agg({'p': np.mean, 'logp': np.mean, 'count': np.sum})
            if len(temp[temp['count'] != num_bg]) > 0:
                print(mirname, ix)
                raise ValueError(f'expected {num_bg} background sequences')
            if len(temp) != len(seqs):
                print(mirname, ix, len(temp), len(seqs))
                raise ValueError('not all seqs')

            SA_bg.append(temp.drop('count', 1))
        SA_bg = pd.concat(SA_bg).reset_index()
        SA_bg_X = [SA_bg]

        # add 12mer sequences for edge sequences
        for ix in range(3):
            temp1 = SA_bg.copy()
            temp1['12mer'] = [('X' * (ix + 1)) + x[ix + 1:] for x in temp1['12mer']]
            temp1 = temp1.groupby(['12mer']).agg(np.mean).reset_index()

            temp2 = SA_bg.copy()
            temp2['12mer'] = [x[:-(ix + 1)] + ('X' * (ix + 1)) for x in temp2['12mer']]
            temp2 = temp2.groupby(['12mer']).agg(np.mean).reset_index()
            SA_bg_X += [temp1, temp2]

        SA_bg_X = pd.concat(SA_bg_X)
        SA_bg_X['mir'] = mirname
        SA_bg_X['stype'] = [utils.get_centered_stype(site8, x) for x in SA_bg_X['12mer'].values]
        print(mirname, len(SA_bg_X[SA_bg_X['stype'] != 'no site']))
        SA_bg_X.to_csv(os.path.join(outfile.replace('MIR', mirname)), sep='\t', index=False, float_format='%.4f')

def combine_res(MIRNAME, MIRSEQ, INFILE_SEQS, INFILE_BG, OUTFILE, MIRDATA=None, NBINS=10, NUM_BG=200, PASSENGER=None):

    if MIRSEQ is not None: #Added to simplify code to one sequence with no passenger
        mirnames = [MIRNAME]
        site8s = [utils.rev_comp(MIRSEQ[1:8]) + 'A']

        # read in RNAplfold results
        iter_res(mirnames, site8s, NBINS, INFILE_SEQS, INFILE_BG, NUM_BG, OUTFILE)

    else:
        mirseqs = pd.read_csv(MIRDATA, sep='\t', index_col='mir')

        for row in mirseqs.iterrows():
            # get names and 8mer sites for the guide and passenger strands, if applicable
            mirnames = [row[0]]
            site8s = [utils.rev_comp(row[1]['guide_seq'][1:8]) + 'A']
            if PASSENGER:
                mirnames += [row[0] + '_pass']
                site8s += [utils.rev_comp(row[1]['pass_seq'][1:8]) + 'A']

            # read in RNAplfold results
            iter_res(mirnames, site8s, NBINS, INFILE_SEQS, INFILE_BG, NUM_BG, OUTFILE)
