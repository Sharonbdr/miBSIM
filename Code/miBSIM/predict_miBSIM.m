function [] = predict_miBSIM(mirname, mirseq, transcript_file, job_name)

if ~exist('../Output', 'dir')
   mkdir('../Output');
end

% Data files
genes = readtable(['../',transcript_file]);
genes = format_gene_table(genes);
miRs = containers.Map(mirname,mirseq);
bio_results_path='../Features/Biochemical/reports/bio_output/';

%% Find potential binding sites and arrange in table
Input = table(0);
keys = miRs.keys;
keep_all=0; %keep_all- whether or not we should keep all possible sites, or include only 3'UTR sites and sites that reside in the last third of the orf

cd Seeds/
for i = 1:height(genes)
    UTR5 = genes.UTR5{i};
    ORF = genes.ORF{i};
    UTR3 = genes.UTR3{i};
    RNA = [UTR5, ORF, UTR3];
    UTR5_end = length(UTR5);
    ORF_end = length(UTR5) + length(ORF);
    
    for j = 1:length(keys)
        miR_name = (keys{j});
        miR = miRs(miR_name);
        temp_table = find_potential_seeds(RNA,miR,UTR5_end,ORF_end,'canonical',keep_all); %CHECKED
        if height(temp_table) > 0
            temp_table.gene_num = i * ones(height(temp_table),1);
            temp_table.miRNA = cellstr(repmat(miR_name,height(temp_table),1));
            if width(Input) == 1
                Input = temp_table;
            else
                Input = [Input; temp_table]; 
            end
        end
    end
end

save(sprintf('../../Output/Input_%s.mat',job_name),'Input','-v7.3');

%% Get Nominal miRNA-mediated repression (no binding site interaction)
cd ../Integrative/

% Calc base features
base_output = base_feat_table(miRs,genes,Input,job_name); 

save(sprintf('../../Output/Base_features_%s.mat',job_name),'base_output','-v7.3');

% Combine all features
int_features = add_bio_features(genes,base_output,cnn_results_path); 

% predict Integrative model
int_pred = predict_Intmodel(int_features); 

save(sprintf('../../Output/Integrative_output_%s.mat',job_name),'features');

%% Get miBSIM repression
cd ../miBSIM/

% Adds neighboring site features
int_pred = prepData(int_pred); 

save(sprintf('../../Output/Neighboring_int_pred_%s.mat',job_name),'int_pred');

% Predicts repression by miBSIM
% **NOTE** to get full feature set and per binding site repression return
% second output and save it.
[predictions,~] = correct_pred(int_pred,genes); 

save(sprintf('../../Output/%s_miBSIM_pred.mat',job_name),'predictions');


end



function genes_in = format_gene_table(genes_in)
for i=1:height(genes_in)
    genes_in.phastcons100{i} = [str2num(genes_in.phastcons100{i})]';
end
for i=1:height(genes_in)
    genes_in.phastcons20{i} = [str2num(genes_in.phastcons20{i})]';
end
for i=1:height(genes_in)
    genes_in.phylops100{i} = [str2num(genes_in.phylops100{i})]';
end
for i=1:height(genes_in)
    genes_in.phylops20{i} = [str2num(genes_in.phylops20{i})]';
end
end