#!/bin/bash


# ==============================================================
# Creating transcript features for the biochemical model
# ==============================================================


# Local Paths
# ==============================================================

# miRNA_models Model
mirna_models_path="${ROOT_path}Biochemical/NN/" 

# Code path and files
explore_path="${ROOT_path}Biochemical"
report_path="${explore_path}/reports" 
input_path="${report_path}/inputs/${dataset}" 
output_path="${report_path}/outputs" 
orf_utr3_fasta="${input_path}/orf_utr3.fa" 
transcripts_file="${input_path}/transcripts.txt" 
rnafld_path="${output_path}/rnaplfold/rnaplfold_orf_utr3/" 

# ==============================================================
# RNAplfold ORF and 3'UTR
# ==============================================================

echo "RNAplfold ORF and UTR3 in $orf_utr3_fasta fasta file..."

mkdir "/scratch/$user/" 
cd "/scratch/$user/" 
RNAplfold -L 40 -W 80 -u 15 < $orf_utr3_fasta

echo "Done."

# ==============================================================
# Processing results
# ==============================================================


echo "Processing results..."

cd $mirna_models_path

python rnaplfold/process_mRNA_folding.py \
--transcripts $transcripts_file \
--indir /scratch/$user/ \
--outdir $rnafld_path \

rm -r "/scratch/$user/" 
echo "Done."


cd $explore_path