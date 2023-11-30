#!/bin/bash

# Paths and Variables
mirname="mir138"
mirseq="AGCTGGTGTTGTGAATCAGGCCG"
transcript_file="Data/transcripts.txt"
job_name="ex"
py_env="/miniconda3/envs/miBSIM_env" 


# Get Biochemical features
source activate $py_env
python Code/Features/Biochemical/get_biofeats.py \
--name ${mirname} \
--mirseq ${mirseq} \
--transcripts ${transcript_file} \
--job_name ${job_name} \                                 
--remove_files True \


# miBSIM model
cd Code/
matlab -nojvm -nodesktop -nodisplay -singleCompThread -r "predict_miBSIM('$mirname','$mirseq','$transcript_file','$job_name'); quit;"


echo DONE
