function [cons_binding,cons_binding_seed,cons_duplex,cons_duplex_seed] = ...
    calc_cons_thermo(UTR5,ORF,UTR3,miRNA,miR_start,RNA_pos1,RNA_start,seed_type,seed_length,glob_path,vienna_dir)
%Calculating the mean dg_binding and dg_duplex of randomized ORF

%RNA = The whole RNA strand
%miRNA = The whole miRNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length of seed

addpath('../../Utils')
perm_num = 1000;


cd('../Thermo')
cons_binding = zeros(1,perm_num);
cons_binding_seed = zeros(1,perm_num);
cons_duplex = zeros(1,perm_num);
cons_duplex_seed = zeros(1,perm_num);

for i = 1:perm_num
    cd('../../Utils')
    random_ORF = permute_orf(ORF);
    cd('../Features/Thermo')
    random_RNA = [UTR5,random_ORF,UTR3];
    [cons_duplex(i),cons_duplex_seed(i),cons_binding(i),cons_binding_seed(i)] =...
        calc_dg_duplex_binding(random_RNA,miRNA,miR_start,RNA_pos1,RNA_start,seed_type,seed_length,glob_path,vienna_dir);
end

cd('../Conservation')

cons_binding = mean(cons_binding);
cons_binding_seed = mean(cons_binding_seed);
cons_duplex = mean(cons_duplex);
cons_duplex_seed = mean(cons_duplex_seed);