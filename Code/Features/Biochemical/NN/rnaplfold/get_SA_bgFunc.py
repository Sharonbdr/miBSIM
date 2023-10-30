from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import shutil
import subprocess
import time
import platform

import numpy as np
import pandas as pd

import utils


def get_SA_bg(SEQUENCE_FILE, TEMP_FOLDER, OUTFILE, NUM_BG=200, NUM_PROCESSES=1, VIEPATH=None):

    if (not os.path.isdir(TEMP_FOLDER)):
        os.makedirs(TEMP_FOLDER)

    SEQUENCES = pd.read_csv(SEQUENCE_FILE, sep='\t', chunksize=NUM_PROCESSES)

    def get_RNAplfold(seq):
        cwd = os.getcwd()
        temp_subfolder = os.path.join(TEMP_FOLDER, seq)
        if os.path.isdir(temp_subfolder):
            shutil.rmtree(temp_subfolder)
        os.makedirs(temp_subfolder)
        os.chdir(temp_subfolder)

        # write sequence to a temporary file
        with open('temp.fa', 'w') as f:
            for ix in range(NUM_BG):
                bg = utils.generate_random_seq(14) + seq + utils.generate_random_seq(14)
                f.write('>{}\n{}\n'.format(ix, bg))

        # # call RNAplfold
        t0 = time.time()
        if VIEPATH is None:
            if platform.system()=="Windows":
                error("Windows system needs a direct path to ViennaRNA installation")
            else:
                mycall = ['RNAplfold', '-L', '40', '-W', '40', '-u', '12']
        else:
            mycall = [VIEPATH+'RNAplfold', '-L', '40', '-W', '40', '-u', '12']

        with open('temp.fa', 'r') as f:
            subprocess.call(mycall, shell=False, stdin=f, stdout=subprocess.PIPE)

        bg_vals = []
        for ix in range(NUM_BG):
            lunp_file = '{}_lunp'.format(ix)
            try:
                rnaplfold_data = pd.read_csv(lunp_file, sep='\t', header=1).set_index(' #i$').astype(float)
            except:
                print('System does not process ViennaRNA properly, check installation or correct path to program.')
            bg_vals.append(rnaplfold_data.loc[26]['12'])

        os.chdir(cwd)
        shutil.rmtree(temp_subfolder)

        return bg_vals

    total = 0
    with open(OUTFILE, 'w') as outfile:
        for chunk in SEQUENCES:
            T0 = time.time()
            seqs = list(chunk['12mer'].values)
            total += len(seqs)
            print('Processing {}'.format(total))

            #with ProcessPoolExecutor(max_workers=options.NUM_PROCESSES) as executor:
                #results = executor.map(get_RNAplfold, seqs)

                #for seq, result in zip(seqs, results):
                    #for res in result:
                        #outfile.write('{}\t{:.7f}\t{:.4f}\n'.format(seq, res, np.log(res)))

            results = [get_RNAplfold(seqs[0])]

            for seq, result in zip(seqs, results):
                for res in result:
                    outfile.write('{}\t{:.7f}\t{:.4f}\n'.format(seq, res, np.log(res)))

            print(time.time() - T0)


