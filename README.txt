This repository contains models for predicting miRNA-mediated repression predictedd by our miBSIM tool. 
system requirments are listed in requirments.txt

## Input Data
This tool requires a folder in Data containing a genes.mat file and miRs.mat file:

	genes :  a table of mRNAs containing the following columns- Name, UTR5, ORF, UTR3, phastcons100, phylops100, phastcons20, phylops20. 
		   Conservation vectors can be obtained using the UCSC Genome Browser Database [https://hgdownload.soe.ucsc.edu/downloads.html under Human genomes, hg38, Multiple alignment].
	miRS  :  a containres.Map object containing miRNA names as keys, and their sequence as values. 

An input example is provided in the Data/ex folder.


## Script
The final model utilizes several published tools, therefore it requires numerous steps to calculate the final repression prediction. 

1. find potential binding sites and create general feature list (excluding biochemical features).
	This is done using the MATLAB functions:
		- Code/Integrative/create_Input('dataset')
		- Code/Integrative/base_feat_table('dataset','user')

2. Calculate biochemical features.
	Biochemical features are calculated using tools published by McGeary et al [https://github.com/kslin/miRNA_models].
	The following scripts should be run:
		- Code/Features/Biochemical/get_transcript_Biofeats.sh
		- Code/Features/Biochemical/get_miR_Biofeats.sh
		- Code/Features/Biochemical/predict_biochemical_feats.sh 

* NOTE: miRNA biochemical script is advised to be run in parallel or completely skipped as it requires a high run time. Please refer to Data/existing_mir_list.txt to see if miRNA was pre-processed and does not require script rerun. 

3. Integrative model prediction.
	This is done using the MATLAB function:
		- Code/Integrative/get_nominal_repression('dataset')

4. miBSIM prediction.
	Corrects former per-site repression prediction to include neighboring site information, by calling the MATLAB function:
		- Code/miBSIM/get_miBSIM_repression('dataset')

An example bash script is provided under example_script.sh