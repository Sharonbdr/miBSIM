function Output = base_feat_table(miRs,genes,Input,user)
% Calculates Base features set (all Integrative model features excluding
% Biochemical features).
%
% dataset    =  A folder in Data that contains a genes table and miRs
%               Map Container.
% user       =  user string, for the purpose of creating a scratch dir.
%


cd ../Features/

Output = Input;
for i = 1:height(Output)
    gene_num = Output.gene_num(i);
    miRNA = upper(miRs(Output.miRNA{i}));
    miRNA(miRNA == 'T') = 'U';
    
    curr_feats = get_feats(genes.UTR5{gene_num},genes.ORF{gene_num},genes.UTR3{gene_num},Output.region{i},Output.seed_length(i),...
        Output.RNA_start(i),miRNA,Output.seed_type{i},user,1,genes.phastcons100{gene_num},...
        genes.phastcons20{gene_num},genes.phylops100{gene_num},genes.phylops20{gene_num});
    
    if i == 1
        all_feats = curr_feats;
    else
        all_feats = [all_feats; curr_feats];
    end
end

Output = [Output, all_feats];


end
