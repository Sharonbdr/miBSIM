function [phastcons_seed100,phastcons_site100,phastcons_flank100,phylops_seed100,phylops_site100,phylops_flank100,...
    phastcons_seed20,phastcons_site20,phastcons_flank20,phylops_seed20,phylops_site20,phylops_flank20] =...
    calc_cons(phast100_vec,phast20_vec,phylop100_vec,phylop20_vec,RNA_start,seed_length)
%Extracting phastCons/phyloP values from the gene's phastCons/phyloP vector
%Parameters:
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7

seed_coor = RNA_start + seed_length-1 - 5:RNA_start + seed_length - 1;
site_coor = RNA_start:RNA_start + seed_length - 1;
flank_start = max(RNA_start-40,1);
flank_end = min(RNA_start+seed_length+39,length(phast100_vec));
flank_coor = [flank_start:RNA_start-1,RNA_start+seed_length:flank_end];

phastcons_seed100 = mean(phast100_vec(seed_coor));
phastcons_site100 = mean(phast100_vec(site_coor));
phastcons_flank100 = mean(phast100_vec(flank_coor));

phastcons_seed20 = mean(phast20_vec(seed_coor));
phastcons_site20 = mean(phast20_vec(site_coor));
phastcons_flank20 = mean(phast20_vec(flank_coor));

phylops_seed100 = phylop100_vec(seed_coor);
phylops_site100 = mean(phylop100_vec(site_coor));
phylops_flank100 = mean(phylop100_vec(flank_coor));

phylops_seed20 = phylop20_vec(seed_coor);
phylops_site20 = mean(phylop20_vec(site_coor));
phylops_flank20 = mean(phylop20_vec(flank_coor));