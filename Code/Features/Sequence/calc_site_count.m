function [utr3_6mer_A1, orf_6mer_A1, utr3_offset_7mer, orf_offset_7mer, utr3_offset_6mer, orf_offset_6mer, ...
    utr3_centered, orf_centered] = calc_site_count(RNA,miR,UTR5_end,ORF_end)

addpath('../../Seeds')

seed_table = find_potential_seeds(RNA,miR,UTR5_end,ORF_end,'all');
if isempty(seed_table)
    utr3_6mer_A1 = 0;
    orf_6mer_A1 = 0;
    utr3_offset_7mer = 0;
    orf_offset_7mer = 0;
    utr3_offset_6mer = 0;
    orf_offset_6mer = 0;
    utr3_centered = 0;
    orf_centered = 0;
    else
    utr3_6mer_A1 = sum(strcmp(seed_table.seed_type,'6mer-A1') & strcmp(seed_table.region,'UTR3'));
    orf_6mer_A1 = sum(strcmp(seed_table.seed_type,'6mer-A1') & strcmp(seed_table.region,'ORF'));
    utr3_offset_7mer = sum(strcmp(seed_table.seed_type,'offset-7mer') & strcmp(seed_table.region,'UTR3'));
    orf_offset_7mer = sum(strcmp(seed_table.seed_type,'offset-7mer') & strcmp(seed_table.region,'ORF'));
    utr3_offset_6mer = sum(strcmp(seed_table.seed_type,'offset-6mer') & strcmp(seed_table.region,'UTR3'));
    orf_offset_6mer = sum(strcmp(seed_table.seed_type,'offset-6mer') & strcmp(seed_table.region,'ORF'));
    utr3_centered = sum(strcmp(seed_table.seed_type,'centered') & strcmp(seed_table.region,'UTR3'));
    orf_centered = sum(strcmp(seed_table.seed_type,'centered') & strcmp(seed_table.region,'ORF'));
end

end