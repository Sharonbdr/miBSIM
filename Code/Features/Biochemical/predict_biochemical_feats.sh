#!/bin/bash

# ===============================================================
# Running Biochemical+ model to get biochemical features
# ===============================================================

# Local Paths
# ===============================================================

# miRNA_models Model
mirna_models_path="${ROOT_path}Biochemical/NN/" 

# Code path and files
explore_path="${ROOT_path}Biochemical" 
report_path="${explore_path}/reports" 
input_path="${report_path}/inputs/${dataset}" 
output_path="${report_path}/outputs"
SA_path="${report_path}/outputs/SA_background" 
mirseq_file="${input_path}/mirseqs.txt"
bg_vals_processed_path="${SA_path}/bg_vals_processed"  
kds_path="${report_path}/outputs/kds" 
features_path="${output_path}/features/${dataset}" 
transcripts_file="${input_path}/transcripts.txt" 
rnafld_path="${output_path}/rnaplfold/rnaplfold_orf_utr3/"
mir_list="${input_path}/miRs.txt"
bio_chemical_model="${mirna_models_path}/biochem_model/trained_models/biochemplus.json"
pred_path="${output_path}/predictions/${dataset}"

# Input parameters
# ===============================================================

overlap_dist=12
upstream_limit=15 
freeago=-7
freeago_pass=-20

# ---------------------------------------------------------------

echo "Creating features..."

cd $mirna_models_path

for mirname in $(cat $mir_list) ; do \
echo "processing $mirname..."
python get_features/write_sites.py \
--transcripts $transcripts_file \
--mir "$mirname" \
--mirseqs $mirseq_file \
--kds $kds_path/"$mirname"_kds.txt \
--sa_bg $bg_vals_processed_path/canon_"$mirname"_bg_vals.txt \
--rnaplfold_dir $rnafld_path \
--overlap_dist $overlap_dist \
--upstream_limit $upstream_limit \
--outfile $features_path/"$mirname".txt ; \
done

echo "Done."

echo "Predicting mRNA repression..."

for mi in $(cat $mir_list) ; do \
echo "processing $mi..."
python biochem_model/predict.py \
--features $features_path/"$mi".txt \
--model $bio_chemical_model \
--freeAGO $freeago \
--outfile $pred_path/"$mi".txt ; \
done

echo "Done."

# Removes partial miRNA files
cd $input_path

if  [ -f mirseqs_rel.txt ]
then
  rm mirseqs_rel.txt
fi

if  [ -f miRs_rel.txt ]
then
  rm miRs_rel.txt
fi


cd $explore_path



