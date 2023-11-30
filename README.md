---
Overview
---
This repository contains model applications for predicting miRNA-mediated repression predicted by our miBSIM tool. 

---
System requirements
---

- MATLAB 2021b and up
- python 3.6 and up
- python virtual environment listed in the requirements.txt file


---
Data input
---
This tool requires a txt file with transcripts of interest, columns delimited by '\t', vectors by ',':

 - transcripts :  a table of mRNAs containing the following columns- Name, UTR5, ORF, UTR3, phastcons100, phylops100, phastcons20, phylops20. 
   Conservation vectors can be obtained using the [UCSC Genome Browser Database](https://hgdownload.soe.ucsc.edu/downloads.html) under Human genomes, hg38, Multiple alignment.

An example is provided at Data/genes.txt.

---
Script
---
The final model utilizes several published tools, therefore it requires numerous steps to calculate the final repression prediction. 

1. Get biochemical features from the published Biochemical+ tool: 
	Biochemical features are calculated using tools published by [McGeary et al](https://github.com/kslin/miRNA_models).
	Published code has been modified to allow for ease of execution, and according to the specific parameters used in our paper.
	This is done using the python function:

		python Code/Features/Biochemical/get_biofeats.py \
		--name mirname \
		--mirseq mirseq \
		--transcripts transcript_file \
		--job_name job_name \ 
		--remove_files False \   # optional
		--skip_mir False \       # optional                                

* **NOTE**: since miRNA biochemical feature generation is streamlined and sequential, step 1 computation time may be extensive. 
Therefore, we added an option to skip redundant feature generation in the case these have already been generated for specific miRNA in the past.
In this case, please make sure to save sub-proccesses files under the --remove_files **False** optional flag to prevent the miRNA feature files to be deleted,
(defaulted to **False**) and to skip re-running these calculations by setting --skip_mir to **True**. These files are saved at the bio_output directory.
This is relevant in the case of running different subsets of transcripts for the same miRNA, and allows for a much faster run-time.


2. Get remaining features and predict miBSIM: 
	Generates remaining features, uses Integrative model to predict nominal repression and finaly corrects this prediction to include
	binding site interaction using the miBSIM model. This is done by calling the MATLAB function:

		predict_miBSIM('mirname','mirseq','transcript_file','job_name')

An example bash script is provided under example_script.sh

