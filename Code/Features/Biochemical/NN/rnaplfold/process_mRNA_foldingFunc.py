import os

import numpy as np
import pandas as pd


def mrna_fold(TRANSCRIPTS, INDIR, OUTDIR):

    transcripts = pd.read_csv(TRANSCRIPTS, sep='\t', index_col=0)

    for row in transcripts.iterrows():
        transcript_length = len(row[1]['ORF']) + len(row[1]['UTR3'])

        # read in rnaplfold outputs
        lunp_file = os.path.join(INDIR, row[0]) + '_lunp'
        rnaplfold_data = pd.read_csv(lunp_file, sep='\t', header=1)
        rnaplfold_data = rnaplfold_data[rnaplfold_data.columns[:15]]

        # use parameters from Agarwal et al., 2015
        rnaplfold_data.columns = ['end'] + list(np.arange(14) + 1)
        rnaplfold_data = rnaplfold_data.set_index('end').astype(float)

        for ix in range(13):
            temp = rnaplfold_data.loc[ix+1]
            rnaplfold_data.loc[ix+1] = rnaplfold_data.loc[ix+1].fillna(temp[ix + 1])

        new_row = pd.Series({14: rnaplfold_data.iloc[-1][13]}, name=transcript_length+1)
        rnaplfold_data = rnaplfold_data[[14]].append(new_row)
        assert(len(rnaplfold_data) == len(rnaplfold_data.dropna()))

        # write to outfile
        outfile = os.path.join(OUTDIR, f'{row[0]}.txt')
        rnaplfold_data.to_csv(outfile, sep='\t')
