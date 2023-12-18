import shutil
from optparse import OptionParser
import platform
import os
import sys
import pandas as pd
import subprocess

def main_prep():
    # Create Biochemical Output directory
    od = os.getcwd() + '/bio_output/'
    md = od + 'util_files/'
    if (not os.path.isdir(md)):
        os.makedirs(md)

    # Create scratch folder
    scdir = os.getcwd() +  '/scratch/'
    if (not os.path.isdir(scdir)):
        os.makedirs(scdir)

    return od, md, scdir

def create_bio_input(transcripts_file, od):
    # Create fasta orf_utr3 file
    trans = pd.read_csv(transcripts_file, sep='\t', index_col=0)

    with open(od+'orf_utr3.fa', 'w') as f:
        for index, row in trans.iterrows():
            f.write(f">{index}\n{row['ORF']+row['UTR3']}\n")



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--name", dest="MIRNAME", help="miRNA name")
    parser.add_option("--mirseq", dest="MIRSEQ", help="miRNA sequence, use when only one sequence and no passenger",
                      default=None)
    parser.add_option("--transcripts", dest="TRANSCRIPTS", help=" file path with transcript sequences")
    parser.add_option("--vienna_path", dest="VIEPATH", help="path to ViennaRNA installation")
    parser.add_option("--job_name", dest="JOB", help="unique job name", default="")
    parser.add_option('--remove_files', dest='RMDIR', action='store_true', default=False, help="remove feature sub process files after biochemical+ prediciton")
    parser.add_option('--skip_mir', dest='SKIP', action='store_true', default=False, help="skip miRNA feature generation (assumes miRNA feature files have been already created)")

    (options, args) = parser.parse_args()


    if options.VIEPATH is None:
        if platform.system() == "Windows":
            raise Warning("Windows may need a direct path to ViennaRNA installation, proceeding with algorithm")
            vienna_path = '../Thermo/ViennaRNA/'
        else:
            vienna_path = '../Thermo/ViennaRNA/'
    else:
        vienna_path = options.VIEPATH
    

    if options.TRANSCRIPTS is None:
        options.TRANSCRIPTS = os.getcwd() + '/Data/transcripts.txt'
    else:
        options.TRANSCRIPTS = os.getcwd() + '/' + options.TRANSCRIPTS

    # Create directories and additional input files
    os.chdir('Code/Features/Biochemical/')
    outdir, utldir, temp_folder = main_prep()
    create_bio_input(options.TRANSCRIPTS, utldir)

    sys.path.append(os.getcwd() + '/NN/')
    sys.path.append(os.getcwd()+'/NN/cnn/')
    sys.path.append(os.getcwd() + '/NN/rnaplfold/')
    sys.path.append(os.getcwd() + '/NN/get_features/')
    sys.path.append(os.getcwd() + '/NN/biochem_model/')
    print('Created input and folders')

    ## Process miRNA features
    if not options.SKIP:
        os.chdir('NN/cnn')
        import generate_12mer_kdsFunc as gen12
        gen12.call_kds(MIRNAME = options.MIRNAME, MIRSEQ = options.MIRSEQ, MIRLEN = 10, LOAD_MODEL = "trained_model/model-100",
                 OUTFILE = utldir+options.MIRNAME+'_kds.txt')
        print('kds')
        os.chdir('../rnaplfold')
        import partition_seqsFunc as part_s
        part_s.part_seqs(MIRNAME = options.MIRNAME, MIRSEQ = options.MIRSEQ, OUTDIR = utldir)
        print('partitioned')
        from get_SA_bgFunc import *
        for i in range(10):
            get_SA_bg(SEQUENCE_FILE = utldir + 'canon_'+options.MIRNAME+'_%d.txt'%i, TEMP_FOLDER = temp_folder,
                      OUTFILE = utldir + 'canon_'+options.MIRNAME+'_%d_bg_vals.txt'%i)
            print('done %i out of 9'%i)
        from combine_resultsFunc import *
        combine_res(MIRNAME = options.MIRNAME, MIRSEQ = options.MIRSEQ, INFILE_SEQS = utldir + 'canon_MIR_IX.txt',
                    INFILE_BG = utldir + 'canon_MIR_IX_bg_vals.txt', OUTFILE = utldir+'canon_'+options.MIRNAME+'bg_vals.txt')
        print('done miRNA')


    ## Process transcripts features
    os.chdir(temp_folder)
    shutil.copyfile('../../Thermo/ViennaRNA/RNAplfold', temp_folder + '/RNAplfold')
    os.chmod(temp_folder + '/RNAplfold', 0o755)
    mycall = [temp_folder + '/RNAplfold', '-L', '40', '-W', '80', '-u', '15']
    with open(utldir+'orf_utr3.fa', 'r') as f:
        subprocess.call(mycall, shell=False, stdin=f, stdout=subprocess.PIPE)
    print('done mRNA folding')

    # process
    os.chdir('../NN/rnaplfold')
    from process_mRNA_foldingFunc import *
    mrna_fold(TRANSCRIPTS = options.TRANSCRIPTS, INDIR = temp_folder, OUTDIR = utldir)
    print('done mRNA')

    ## Apply Biochemical+ model
    os.chdir('../get_features')
    from write_sitesFunc import *
    write_sites(MIRNAME = options.MIRNAME, MIRSEQ_IN = options.MIRSEQ, TRANSCRIPTS_IN = options.TRANSCRIPTS,
                KDS_IN = utldir+options.MIRNAME+'_kds.txt', SA_BG = utldir+'canon_'+options.MIRNAME+'bg_vals.txt',
                RNAPLFOLD_DIR = utldir, OUTFILE = utldir+options.JOB+'_'+options.MIRNAME+'_feats.txt')
    print('done sites')

    os.chdir('../biochem_model')
    from predictFunc import *
    predict(FEATURES = utldir+options.JOB+'_'+options.MIRNAME+'_feats.txt', MODEL = 'biochemplus.json',
            OUTFILE = outdir+options.JOB+'_'+options.MIRNAME+'_pred.txt')
    print('finished biochemical processing')

    # Remove processing files (per users request) and temp folder
    if options.RMDIR:
        shutil.rmtree(utldir)
        print("removed files")
    shutil.rmtree(temp_folder)


