function Input = create_Input(dataset)
% Creating an Input table for base_feat_table
%
% dataset    =  A folder in Data that contains a genes table and miRs
%               Map Container.
%

% Variables
load(sprintf('../../Data/%s/genes.mat',dataset),'genes'); 
load(sprintf('../../Data/%s/miRs.mat',dataset),'miRs');

Input = table(0);
keys = miRs.keys;
keep_all=0; %keep_all- whether or not we should keep all possible sites, or include only 3'UTR sites and sites that reside in the last third of the orf

% Creates directories
if ~exist('../../Output', 'dir')
   mkdir('../../Output');
end
if ~exist(sprintf('../Features/Biochemical/reports/inputs/%s',dataset), 'dir')
   mkdir(sprintf('../Features/Biochemical/reports/inputs/%s',dataset));
end
if ~exist('../Features/Biochemical/reports/outputs/rnaplfold', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/rnaplfold');
end
if ~exist('../Features/Biochemical/reports/outputs/rnaplfold/rnaplfold_orf_utr3', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/rnaplfold/rnaplfold_orf_utr3');
end
if ~exist('../Features/Biochemical/reports/outputs/SA_background/bg_vals', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/SA_background/bg_vals');
end
if ~exist('../Features/Biochemical/reports/outputs/SA_background/sequences', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/SA_background/sequences');
end
if ~exist('../Features/Biochemical/reports/outputs/features', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/features');
end
if ~exist(sprintf('../Features/Biochemical/reports/outputs/features/%s',dataset), 'dir')
   mkdir(sprintf('../Features/Biochemical/reports/outputs/features/%s',dataset));
end
if ~exist('../Features/Biochemical/reports/outputs/predictions', 'dir')
   mkdir('../Features/Biochemical/reports/outputs/predictions');
end
if ~exist(sprintf('../Features/Biochemical/reports/outputs/predictions/%s',dataset), 'dir')
   mkdir(sprintf('../Features/Biochemical/reports/outputs/predictions/%s',dataset));
end

% Find potential binding sites and arrange in table
cd ../Seeds/
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
        temp_table = find_potential_seeds(RNA,miR,UTR5_end,ORF_end,'canonical',keep_all);
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


save(sprintf('../../Output/Input_%s.mat',dataset),'Input','-v7.3');


%% Preps data for Biochemical features:
% missing miRNAs are mirnas that were not pre-processed (list under Data/existing_mir_list)
ext_mir=readtable('../../Data/existing_mir_list.txt');
missing_mirs=[];
for i=1:length(keys)
    curr_mir=lower(replace(keys{i},'-','_'));
    ind=find(strcmp(ext_mir.mirname,curr_mir));
    if length(ind)<1
        missing_mirs=[missing_mirs;i];
    end
end


cd ../Features/Biochemical/reports/inputs/
create_BioInput_trans(dataset);
create_BioInput_mirs(dataset,missing_mirs);

cd ../../../../Integrative/


end
