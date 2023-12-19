---
miBSIM - miRNA Binding Site Interaction Model
---
This repository contains model applications for predicting miRNA-mediated repression predicted by our miBSIM tool. 

---
System requirements
---

- MATLAB 2021b and up
- python 3.6
- python virtual environment listed in the requirements.txt file


---
Data input
---
This tool requires a txt file with transcripts of interest, columns delimited by '\t', vectors by ',':

 - transcripts :  a table of mRNAs containing the following columns- Name, UTR5, ORF, UTR3, phastcons100, phylops100, phastcons20, phylops20. 
   Conservation vectors can be obtained using the [UCSC Genome Browser Database](https://hgdownload.soe.ucsc.edu/downloads.html) under Human genomes > hg38 > Multiple alignment.

An example is provided at Data/transcripts.txt.

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
		--remove_files \#     optional
		--skip_mir \#         optional                                

--remove_files	: defaulte = FALSE, option to delete intermediate files created by the biochemical+ model.
--skip_mir	: default = FALSE, option to skip miRNA biochemical feature generation assuming these files exist under the miR_DataFiles directory.

* **NOTE**: miRNA biochemical feature generation was modified from original publication to be streamlined and sequential. This causes step 1 computation time to be extensive. 
Therefore, we added the option to skip miRNA biochemical feature generation by providing the files under the miR_DataFiles directory.
In this case, call the **--skip_mir**  option (defaulted to **False**). This is also relevant in the case of running different subsets of transcripts for the same miRNA, and allows for a much faster run-time.
Preprocessed miRNAs are published in this directory, please make sure to match miRNA name in a case sensitive manner.



2. Get remaining features and predict miBSIM: 
	Generates remaining features, uses Integrative model to predict nominal repression and finaly corrects this prediction to include
	binding site interaction using the miBSIM model. This is done by calling the MATLAB function:

		predict_miBSIM('mirname','mirseq','transcript_file','job_name')

An example bash script is provided under example_script.sh

---
Troubleshooting
---
Please make sure you set all files premissions to read write and execute. This is specificaly important for the following files:

- Code/Features/Thermo/ViennaRNA/*
- Code/Features/Conservation/pmc.markov
- Code/Features/Conservation/pmc.sci
- Code/Features/Conservation/spatt
