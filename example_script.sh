#!/bin/bash

# ========================================================================== #
#                                                                            #
# These models utilize tools publishes by several groups under different     #
# program languages, therefore they require a number of seperate steps.      #
# Code and enviroment requirements are listed in requirements.txt.           #
#                                                                            #
# If running a large dataset it is advised to run some processes in          #
# parallel, namely the Biochemical+ miRNA features.                          #
# Many miRNAs were pre-processes to allow users to skip this step.           #
#                                                                            #
# ========================================================================== #

# ========================================================================== #
#                            Paths and Variables                             #
# ========================================================================== #
cd miBSIM-main
export dataset="ex"
export user="sharon"
export ROOT_path="/tamir2/sharon/miRNA_models/Code/Features/"

python test_py.py
# ========================================================================== #
#                             general features                               #
# ========================================================================== #
# Creates models input files and directories, and calculates general features
#cd Code/Integrative/
#matlab -c 27000@lm8-2 -nojvm -nodesktop -nodisplay -singleCompThread -r "create_Input('$dataset'); base_feat_table('$dataset','$user'); quit;"
#echo "done"

# ========================================================================== #
#                          Biochemical features                              #
# ========================================================================== #
#
# The scripts used here follows: "https://github.com/kslin/miRNA_models"
# published by McGEary et al.
 
# Get Biochemical+ mRNA features
#cd ../Features/Biochemical
#./get_transcript_Biofeats.sh

# Get Biochemical+ miRNA features
# Currently designed to run only if miRNA were not pre-processed
#if  [ -f reports/inputs/$dataset/mirseqs_rel.txt ]
#then
#  ./get_miR_Biofeats.sh
#fi

# predict Biochemical+ features combined
#./predict_biochemical_feats.sh 


# ========================================================================== #
#                             Integrative model                              #
# ========================================================================== #
## Predict Integrative model
#
# combines all features and predicts nominal repression per site (assumes
# sites act independantly)
#

#cd ../../Integrative/
#matlabr2022b -nojvm -nodesktop -nodisplay -singleCompThread -r "get_nominal_repression('$dataset'); quit;"

# ========================================================================== #
#                           predict miBSIM model                             #
# ========================================================================== #
#
# Gets neighboring site information and corrects sites affective repression
#

#cd ../miBSIM/
#matlabr2022b -nojvm -nodesktop -nodisplay -singleCompThread -r "get_miBSIM_repression('$dataset'); quit;"

echo "DONE"
