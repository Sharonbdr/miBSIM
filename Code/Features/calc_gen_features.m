function Output = calc_gen_features(UTR5,ORF,UTR3,region,seed_length,RNA_start,miRNA,seed_type,glob_path,vienna_dir,with_cons_orf,...
    phastcons100,phastcons20,phylops100,phylops20)
% Calculate a 1xn vector of features
%

% UTR5/ORF/UTR3        =   Sequence of the region
% region               =   Of the binding site
% seed_length          =   6/7/8
% RNA_start            =   beginning of the target site in the complete mRNA
% miRNA                =   Sequence of the miRNA
% seed_type            =   '6mer'/'7mer-m8'/'7mer-A1'/'8mer'
% user                 =   user string, for the purpose of creating a scratch dir
% glob_path            =   Path of the folder containing ViennaRNA folder in network
%                          drive
% vienna_dir           =   Path of copied ViennaRNA folder
% with_cons_orf        =   0/1 index for including conservation and ORF features, i.e. if the
%                          gene is human and endogenous or non-human and engineered.
% phastcons and phylop =   Needed if with_cons_orf == 1; phastCons and phyloP score vectors for the gene


addpath('../Utils')
addpath('../Seeds')

UTR5 = upper(UTR5);
UTR5(UTR5 == 'T') = 'U';
ORF = upper(ORF);
ORF(ORF == 'T') = 'U';
UTR3 = upper(UTR3);
UTR3(UTR3 == 'T') = 'U';
RNA = [UTR5, ORF, UTR3];
miRNA = upper(miRNA);
miRNA(miRNA == 'T') = 'U';
ORF_start = length(UTR5) + 1;
ORF_end = length(UTR5) + length(ORF);
if any(strcmp(seed_type,{'6mer','6mer-A1','7mer-A1','7mer-m8','8mer'}))
    miR_start = 2;
else
    miR_start = 3;
end
if strcmp(seed_type,'8mer') || strcmp(seed_type,'7mer-A1')
    RNA_pos1 = RNA_start + seed_length-1;
else
    RNA_pos1 = RNA_start + seed_length;
end

% region_start/end - in local (RNA strand) coordinates
switch region
    case 'UTR5'
        region_start = 1;
        region_end = length(UTR5);
    case 'ORF'
        region_start = length(UTR5)+1;
        region_end = length(UTR5)+length(ORF);
    case 'UTR3'
        region_start = length(UTR5)+length(ORF)+1;
        region_end = length(RNA);
end

Output = table(0);

% Conservation
cd 'Conservation';

Output.prob_binom = calc_prob_binom(RNA,RNA_start,seed_length);
Output.prob_exact = calc_prob_exact(RNA,RNA_start,seed_length);

if with_cons_orf == 1
    [Output.phastcons_seed100,Output.phastcons_site100,Output.phastcons_flank100,phylops_seed100,Output.phylops_site100...
        ,Output.phylops_flank100,Output.phastcons_seed20,Output.phastcons_site20,Output.phastcons_flank20,phylops_seed20,...
        Output.phylops_site20,Output.phylops_flank20] = calc_cons(phastcons100,phastcons20,phylops100,phylops20,RNA_start,seed_length);

    phylops_seed100 = flipud(phylops_seed100);
    for g = 1:length(phylops_seed100)
        Output{1,sprintf('phylops100_pos%d',g)} = phylops_seed100(g);
    end

    phylops_seed20 = flipud(phylops_seed20);
    for g = 1:length(phylops_seed20)
        Output{1,sprintf('phylops20_pos%d',g)} = phylops_seed20(g);
    end
    
end

% Sequence
cd '../Sequence';

[Output.distance_score, Output.dist_start, Output.dist_end] = calc_utr_pos(RNA_start,region_start,region_end,miRNA,seed_length);
Output.au_content = calc_au_content(RNA,RNA_start,seed_length,seed_type);
Output.miR_pos1_A = double(miRNA(1) == 'A');
Output.miR_pos1_C = double(miRNA(1) == 'C');
Output.miR_pos1_G = double(miRNA(1) == 'G');
Output.miR_pos8_A = double(miRNA(8) == 'A');
Output.miR_pos8_C = double(miRNA(8) == 'C');
Output.miR_pos8_G = double(miRNA(8) == 'G');

if RNA_pos1 <= length(RNA)
    Output.RNA_pos1_A = double(RNA(RNA_pos1) == 'A');
    Output.RNA_pos1_C = double(RNA(RNA_pos1) == 'C');
    Output.RNA_pos1_G = double(RNA(RNA_pos1) == 'G');
end

RNA_pos8 = RNA_pos1 - 7;
if RNA_pos8 >= 1
    Output.RNA_pos8_A = double(RNA(RNA_pos8) == 'A');
    Output.RNA_pos8_C = double(RNA(RNA_pos8) == 'C');
    Output.RNA_pos8_G = double(RNA(RNA_pos8) == 'G');
    if RNA_pos8 >= 2
        Output.RNA_pos9_A = double(RNA(RNA_pos8-1) == 'A');
        Output.RNA_pos9_C = double(RNA(RNA_pos8-1) == 'C');
        Output.RNA_pos9_G = double(RNA(RNA_pos8-1) == 'G');
    else
        Output.RNA_pos9_A = 0;
        Output.RNA_pos9_C = 0;
        Output.RNA_pos9_G = 0;
    end
else
    Output.RNA_pos8_A = 0;
    Output.RNA_pos8_C = 0;
    Output.RNA_pos8_G = 0;
    Output.RNA_pos9_A = 0;
    Output.RNA_pos9_C = 0;
    Output.RNA_pos9_G = 0;
end

switch region
    case 'UTR5'
        curr_region = UTR5;
    case 'ORF'
        curr_region = ORF;
    case 'UTR3'
        curr_region = UTR3;
end

Output.dist_start_rel = round(10 ^ Output.dist_start) / length(curr_region);
Output.region_A = count(curr_region,'A') / length(curr_region);
Output.region_C = count(curr_region,'C') / length(curr_region);
Output.region_G = count(curr_region,'G') / length(curr_region);

seed_seq = RNA(RNA_start:RNA_start + seed_length-1);
Output.site_A = count(seed_seq,'A') / seed_length;
Output.site_C = count(seed_seq,'C') / seed_length;
Output.site_G = count(seed_seq,'G') / seed_length;

if ~isempty(UTR5)
    Output.UTR5_len = log10(length(UTR5));
    Output.UTR5_AU = 1 - count(UTR5,{'G','C'}) / length(UTR5);
else
    Output.UTR5_len = 0;
    Output.UTR5_AU = 0;
end
Output.ORF_len = log10(length(ORF));
Output.ORF_AU = 1 - count(ORF,{'G','C'}) / length(ORF);

if ~isempty(UTR3)
    Output.UTR3_len = log10(length(UTR3));
    Output.UTR3_AU = 1 - count(UTR3,{'G','C'}) / length(UTR3);
else
    Output.UTR3_len = 0;
    Output.UTR3_AU = 0;
end

load('../../../Data/miBSIM_param/TS7.mat','TS7');
ind = find(strcmp(TS7.seed,miRNA(2:8)));
Output.TA = TS7.TA(ind);
Output.TA_hela = TS7.TA_hela(ind);
Output.SPS_6mer = TS7.SPS_6mer(ind);
Output.SPS_7mer_m8 = TS7.SPS_7mer_m8(ind);
Output.mean_SPS = TS7.mean_SPS(ind);

miR_end = max(1,RNA_pos1 - (length(miRNA) - 1));
pairing = erase(num2str(is_watson_crick(RNA(miR_end:miR_end+length(miRNA)-1),miRNA)),' ');
tokens = regexp(pairing,'1+','match');
if ~isempty(tokens)
    Output.longest_pairing = max(cellfun(@length,tokens));
else
    Output.longest_pairing = 0;
end

Output.paired_miR_end = sum(is_watson_crick(RNA(miR_end:miR_end+6),miRNA(end-6:end)));
Output.pairing3p_sw = calc_3p_pairing_sw(RNA,RNA_start,miRNA,seed_type);
Output.pairing3p = calc_3p_pairing(RNA,RNA_start,miRNA,seed_length,miR_start);

miR_cog = RNA(miR_end:RNA_pos1);
Output.miR_cog_A = count(miR_cog,'A') / length(miR_cog);
Output.miR_cog_C = count(miR_cog,'C') / length(miR_cog);
Output.miR_cog_G = count(miR_cog,'G') / length(miR_cog);

miR_cog_up = RNA(max(1,miR_end-50):max(1,miR_end-1));
Output.miR_cog_up_A = count(miR_cog_up,'A') / length(miR_cog_up);
Output.miR_cog_up_C = count(miR_cog_up,'C') / length(miR_cog_up);
Output.miR_cog_up_G = count(miR_cog_up,'G') / length(miR_cog_up);
miR_cog_down = RNA(min(length(RNA),RNA_pos1+1):min(length(RNA),RNA_pos1+50));
Output.miR_cog_down_A = count(miR_cog_down,'A') / length(miR_cog_down);
Output.miR_cog_down_C = count(miR_cog_down,'C') / length(miR_cog_down);
Output.miR_cog_down_G = count(miR_cog_down,'G') / length(miR_cog_down);

[ORF_pum,UTR3_pum,nei_pum] = calc_pum(RNA,RNA_start,seed_length,ORF_start,ORF_end);
Output.ORF_pum1 = ORF_pum(1);
Output.ORF_pum2 = ORF_pum(2);
Output.UTR3_pum1 = UTR3_pum(1);
Output.UTR3_pum2 = UTR3_pum(2);
Output.nei_pum1 = nei_pum(1); 
Output.nei_pum2 = nei_pum(2);

RBP_sites = {'UAUUU','UAUUUAU','AUUUAU','GUGAAG','AUUUAUU','UAUUUA','UAUUUU','AGCCA','AGAGAA','UUUAU','UUUUU','UUUUAAA','UUUUUU',...
'UUUUUUU','UUUAAA','UUUAAAA','UUUUG','UGUUU','GAAAA','CAAAA'};
for j = 1:length(RBP_sites)
    Output{1,sprintf('RBP%d',j)} = count(UTR3,RBP_sites{j});
end

seed_table = find_potential_seeds(RNA,miRNA,length(UTR5),ORF_end,'all',1);
seed_table = seed_table(strcmp(seed_table.region,'UTR3'),:);
for seed = ["6mer-A1","offset-7mer","offset-6mer","centered","CDNST1","CDNST2","CDNST3","CDNST4"]
    Output{1,sprintf('utr3_%s',strrep(seed,'-','_'))} = sum(strcmp(seed_table.seed_type,seed));
end

if ~isempty(UTR5)
    UTR5_T = UTR5;
    UTR5_T(UTR5_T == 'U') = 'T';
    Output.utr5_3mer = length(regexp(UTR5_T,seqrcomplement(miRNA(2:4))));
else
    Output.utr5_3mer = 0;
end
if ~isempty(ORF)
    ORF_T = ORF;
    ORF_T(ORF_T == 'U') = 'T';
    Output.orf_3mer = length(regexp(ORF_T,seqrcomplement(miRNA(2:4))));
else
    Output.orf_3mer = 0;
end
if ~isempty(UTR3)
    UTR3_T = UTR3;
    UTR3_T(UTR3_T == 'U') = 'T';
    Output.utr3_3mer = length(regexp(UTR3_T,seqrcomplement(miRNA(2:4))));
else
    Output.utr3_3mer = 0;
end

% Thermodynamic
cd '../Thermo';
[Output.dg_duplex,Output.dg_duplex_seed,Output.dg_binding,Output.dg_binding_seed] = ...
    calc_dg_duplex_binding(RNA,miRNA,miR_start,RNA_pos1,RNA_start,seed_type,seed_length,glob_path,vienna_dir); 

Output.dg_open = calc_dg_open(RNA,RNA_pos1,length(miRNA),glob_path,vienna_dir);

Output.SPS = calc_sps(RNA,miRNA,RNA_start,seed_type,seed_length,glob_path,vienna_dir);

sa = calc_sa(RNA,RNA_start,seed_type,region_end,region,glob_path,vienna_dir);
if ~isnan(sa)
    Output.SA = sa;
else
    Output.SA = 999;
end

max_energy = 0:-5:-30;
mismatch = 2;

MREs = calc_MRE(ORF,miRNA,max_energy,mismatch,glob_path,vienna_dir);
for p = 1:length(MREs)
    Output{1,sprintf('MRE_therm_ORF_m%d_%d',mismatch,abs(max_energy(p)))} = MREs(p);
end

MREs = calc_MRE(UTR3,miRNA,max_energy,mismatch,glob_path,vienna_dir);
for p = 1:length(MREs)
    Output{1,sprintf('MRE_therm_UTR_m%d_%d',mismatch,abs(max_energy(p)))} = MREs(p);
end

% Codons
if with_cons_orf == 1
    cd '../Codons';
    [Output.cub_global_win,Output.cub_global_orf,Output.cub_local_win,Output.cub_local_orf,Output.cub_local_CAI_win...
   ,Output.cub_local_CAI_orf,Output.tAI_score_win,Output.tAI_score_orf,Output.charge_score_win...
   ,Output.charge_score_orf,Output.slow_score_win,Output.slow_score_orf,Output.TDR_score_win,Output.TDR_score_orf] = ...
   calc_codon(RNA,RNA_start,seed_length,ORF_start,ORF_end,region);
end

cd '..'

Output.Var1 = [];

end