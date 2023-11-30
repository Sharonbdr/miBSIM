import os

import numpy as np
import pandas as pd

import get_site_features
import utils

np.set_printoptions(threshold=np.inf, linewidth=200)
pd.options.mode.chained_assignment = None


def write_sites(MIRNAME, MIRSEQ_IN, TRANSCRIPTS_IN, KDS_IN, SA_BG, OUTFILE, RNAPLFOLD_DIR=None, PCT_FILE=None, KD_CUTOFF=np.inf, OVERLAP_DIST=12, UPSTREAM_LIMIT=15, ONLY_CANON=False, MIRDATA=None):

    TRANSCRIPTS = pd.read_csv(TRANSCRIPTS_IN, sep='\t', index_col=0)
    TRANSCRIPTS['orf_length'] = TRANSCRIPTS['ORF'].apply(len)
    TRANSCRIPTS['utr3_length'] = TRANSCRIPTS['UTR3'].apply(len)

    if MIRSEQ_IN is not None:  # Added to simplify code to one sequence with no passenger
        MIRSEQ = MIRSEQ_IN
        FAMILY = MIRNAME
    else:
        mirseqs = pd.read_csv(MIRDATA, sep='\t', index_col='mir')

        if '_pass' in MIRNAME:
            MIRSEQ = mirseqs.loc[MIRNAME.replace('_pass', '')]['pass_seq']
            FAMILY = mirseqs.loc[MIRNAME.replace('_pass', '')]['pass_family']
        else:
            MIRSEQ = mirseqs.loc[MIRNAME]['guide_seq']
            FAMILY = mirseqs.loc[MIRNAME]['guide_family']

    SITE8 = utils.rev_comp(MIRSEQ[1:8]) + 'A'
    print(MIRNAME, SITE8)

    # if KD file provided, find sites based on KD file
    if KDS_IN is not None:
        KDS = pd.read_csv(KDS_IN, sep='\t')

        KDS['mir'] = KDS['mir'].str.lower()
        KDS['mir'] = KDS['mir'].replace('-', '_', regex=True)

        if ONLY_CANON:
            KDS = KDS[KDS['aligned_stype'] != 'no site']
        KDS = KDS[KDS['best_stype'] == KDS['aligned_stype']]

        temp = KDS[KDS['mir'] == FAMILY]
        if len(temp) == 0:
            raise ValueError('{} not in kd files'.format(FAMILY))
        mir_kd_dict = {x: y for (x, y) in zip(temp['12mer'], temp['log_kd']) if (y < KD_CUTOFF)}

        # find all the sites and KDs
        all_features = []
        for row in TRANSCRIPTS.iterrows():
            all_features.append(
                get_site_features.get_sites_from_kd_dict_improved(row[0], row[1]['ORF']+row[1]['UTR3'], mir_kd_dict,
                                                                  OVERLAP_DIST))

    # otherwise, go by sequence
    else:
        all_features = []
        for row in TRANSCRIPTS.iterrows():
            all_features.append(get_site_features.get_sites_from_sequence(row[0], row[1]['ORF']+row[1]['UTR3'], SITE8,
                                                                          overlap_dist=OVERLAP_DIST,
                                                                          only_canon=ONLY_CANON))

    all_features = pd.concat(all_features).sort_values('transcript')
    all_features['mir'] = MIRNAME.replace('_pass', '*')

    # add site accessibility background information
    temp = pd.read_csv(SA_BG, sep='\t', index_col='12mer').reindex(all_features['12mer'].values)

    temp['mir'] = temp['mir'].str.lower()
    temp['mir'] = temp['mir'].replace('-', '_', regex=True)

    all_features['logSA_bg'] = temp['logp'].values

    # add stypes
    all_features['stype'] = [utils.get_centered_stype(SITE8, seq) for seq in all_features['12mer'].values]

    # sanity check on background
    temp = all_features[all_features['stype'] != 'no site']
    if len(temp) != len(temp.dropna()):
        raise ValueError('Error in site accessibility background assignment')

    print('Adding 3p score and SA')

    # add transcript-specific information
    temp = []
    for transcript, group in all_features.groupby('transcript'):
        locs = group['loc'].values

        # add threep pairing score
        sequence = TRANSCRIPTS.loc[transcript]['ORF']+TRANSCRIPTS.loc[transcript]['UTR3']
        group['Threep'] = [
            get_site_features.calculate_threep_score(MIRSEQ, sequence, int(loc - 3), UPSTREAM_LIMIT) for loc in
            locs]

        # add site accessibility information
        lunp_file = os.path.join(RNAPLFOLD_DIR, transcript) + '.txt'
        rnaplfold_data = pd.read_csv(lunp_file, sep='\t', index_col='end')
        group['SA'] = rnaplfold_data.reindex(locs + 7)['14'].values.astype(float)  # Agarwal 2015 parameters
        group['logSA'] = np.log(group['SA'])

        temp.append(group)

    all_features = pd.concat(temp)
    all_features['orf_length'] = TRANSCRIPTS.reindex(all_features['transcript'].values)['orf_length'].values
    all_features['utr3_length'] = TRANSCRIPTS.reindex(all_features['transcript'].values)['utr3_length'].values
    all_features['in_ORF'] = all_features['loc'] < (all_features['orf_length'] + 15)
    all_features['logSA_diff'] = all_features['logSA'] - all_features['logSA_bg']
    all_features['utr3_loc'] = all_features['loc'] - all_features['orf_length']
    all_features['passenger'] = ('_pass' in MIRNAME) or ('*' in MIRNAME)

    print('Adding PCT')

    # add PCT information if indicated
    if PCT_FILE is not None:
        pct_df = pd.read_csv(PCT_FILE, sep='\t',
                             usecols=['Gene ID', 'miRNA family', 'Site type', 'Site start', 'PCT'])
        pct_df.columns = ['transcript', 'mir', 'stype', 'loc', 'PCT']
        pct_df = pct_df[pct_df['mir'] == FAMILY]
        if len(pct_df) == 0:
            all_features['PCT'] = 0
            print(f"No PCT information for {FAMILY}")
        else:
            pct_df['offset'] = [1 if x in ['8mer-1a', '7mer-m8'] else 0 for x in pct_df['stype']]
            pct_df['loc'] = pct_df['loc'] + pct_df['offset']
            pct_df = pct_df[pct_df['stype'] != '6mer']
            pct_df = pct_df.set_index(['transcript', 'loc'])

            temp1 = all_features[all_features['in_ORF']]
            temp1['PCT'] = 0
            temp2 = all_features[~all_features['in_ORF']]
            temp2['PCT'] = pct_df.reindex(temp2[['transcript', 'utr3_loc']])['PCT'].values
            temp2['PCT'] = temp2['PCT'].fillna(0.0)
            all_features = pd.concat([temp1, temp2])

    else:
        print(f"No PCT information for {FAMILY}")
        all_features['PCT'] = 0

    all_features = all_features.set_index('transcript').sort_index()

    # write outputs
    all_features.to_csv(OUTFILE, sep='\t')
