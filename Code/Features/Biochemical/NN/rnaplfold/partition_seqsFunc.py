import os

import numpy as np
import pandas as pd

import utils

def iter_mirs(names, sites, canon, nbins, outdir):
    for mirname, site8 in zip(names, sites):
        all_seqs = utils.generate_12mers(site8, canon)
        print(len(all_seqs))
        temp = pd.DataFrame(
            {'12mer': all_seqs, 'aligned_stype': [utils.get_centered_stype(site8, seq) for seq in all_seqs]})

        # separate by canonical or not canonical and partition into files
        with_site = temp[temp['aligned_stype'] != 'no site']
        with_site['bin'] = [x % nbins for x in np.arange(len(with_site))]
        for ix, group in with_site.groupby('bin'):
            group[['12mer', 'aligned_stype']].to_csv(
                os.path.join(outdir, 'canon_{}_{}.txt'.format(mirname, ix)), sep='\t', index=False)

        if not canon:
            no_site = temp[temp['aligned_stype'] == 'no site']
            no_site['bin'] = [x % nbins for x in np.arange(len(no_site))]
            for ix, group in no_site.groupby('bin'):
                group[['12mer', 'aligned_stype']].to_csv(
                    os.path.join(outdir, 'noncanon_{}_{}.txt'.format(mirname, ix)), sep='\t', index=False)




def part_seqs(MIRNAME, MIRSEQ, OUTDIR, NBINS=10, ONLY_CANON=True, MIRDATA=None, PASSENGER=None):

    if (not os.path.isdir(OUTDIR)):
        os.makedirs(OUTDIR)

    if MIRSEQ is not None: #Added to simplify code to one sequence with no passenger
        mirnames = [MIRNAME]
        mirseq=MIRSEQ
        site8s = [utils.rev_comp(mirseq[1:8]) + 'A']

        iter_mirs(mirnames, site8s, ONLY_CANON, NBINS, OUTDIR)


    else:
        MIRSEQS = pd.read_csv(MIRDATA, sep='\t', index_col='mir')

        for row in MIRSEQS.iterrows():
            # get names and 8mer sites for the guide and passenger strands, if applicable
            mirnames = [row[0]]
            site8s = [utils.rev_comp(row[1]['guide_seq'][1:8]) + 'A']
            if PASSENGER:
                mirnames += [row[0] + '_pass']
                site8s += [utils.rev_comp(row[1]['pass_seq'][1:8]) + 'A']

            iter_mirs(mirnames, site8s, ONLY_CANON, NBINS, OUTDIR)


