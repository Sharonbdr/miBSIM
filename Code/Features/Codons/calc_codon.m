function [cub_global_win,cub_global_orf,cub_local_win,cub_local_orf,cub_local_CAI_win,cub_local_CAI_orf,tAI_score_win,tAI_score_orf,...
    charge_score_win,charge_score_orf,slow_score_win,slow_score_orf,TDR_score_win,TDR_score_orf] = ...
    calc_codon(RNA,RNA_start,seed_length,ORF_start,ORF_end,region)
%Calculating the codon speed scores - window around site (if applicable) and whole ORF

%RNA = The whole RNA strand
%RNA_start = Start position of the site in the RNA
%seed_length = Length between 6 and 7
%ORF_start/end - Start/end position of the site's region in the gene
%region - Region of the seed

load('../../../Data/codon_data.mat') %Includes codon frequency data and TDR & tAI weights

window = 25;

RNA = upper(RNA);
RNA(RNA == 'U') = 'T';

if strcmp(region,'ORF')
    types = {'orf','window'};
    %This is the "relative reading frame", meaning the nt offset needed in order to take
    %the right codons (according to the ORF's reading frame) in the site
    start_relative_rf = mod(RNA_start-ORF_start,3);
    RNA_end = RNA_start + seed_length - 1;
    end_relative_rf = mod(RNA_end-ORF_start,3);
else
    types = {'orf'};
    cub_global_win = 999;
    cub_local_win = 999;
    cub_local_CAI_win = 999;
    tAI_score_win = 999;
    TDR_score_win = 999;
    charge_score_win = 999;
    slow_score_win = 999;
end
    
for t = 1:length(types)
    if strcmp(types{t},'orf')
        win_start = ORF_start;
        win_end = ORF_end; 
    else
        win_start = max(ORF_start, RNA_start - 3 * window - start_relative_rf);
        win_end = min(ORF_end, RNA_end + 2 - end_relative_rf + 3 * window);
    end
    
    window_codons = RNA(win_start:win_end);
    if mod(length(window_codons),3) ~= 0
        error('window not divisible by 3');
    end

    window_aas = nt2aa(window_codons,'AlternativeStartCodons',false);

    win_slow_aas = count(window_aas,{'D','P'});
    slow_score = -1 * win_slow_aas / length(window_aas) + (length(window_aas) - win_slow_aas) / win_slow_aas;
    if isinf(slow_score)
        slow_score = 0;
    end

    win_positive_aas = count(window_aas,{'R','K','H'});
    win_negative_aas = count(window_aas,{'E','D'});
    charge_score = (win_positive_aas - win_negative_aas) / (length(window_aas));

    %Splitting into triplets
    window_codons = regexp(window_codons,'.{3}','match');

    win_scores_global = zeros(length(window_codons),1);
    win_scores_local = zeros(length(window_codons),1);
    win_scores_local_CAI = zeros(length(window_codons),1);
    win_tAI_scores = zeros(length(window_codons),1);
    win_TDR_scores = zeros(length(window_codons),1);

    for c = 1:length(window_codons)
        curr_codon = window_codons{c};
        win_scores_global(c) = CUT_global(curr_codon);
        win_scores_local(c) = CUT_local(curr_codon);
        win_scores_local_CAI(c) = CUT_local_CAI(curr_codon);
        win_tAI_scores(c) = tAI_weights(curr_codon);
        win_TDR_scores(c) = TDR(curr_codon);
    end

    cub_global = geomean(win_scores_global);
    cub_local = geomean(win_scores_local);
    cub_local_CAI = geomean(win_scores_local_CAI);
    tAI_score = geomean(win_tAI_scores);
    TDR_score = geomean(win_TDR_scores(~isnan(win_TDR_scores)));
    if isnan(TDR_score)
        TDR_score = 0;
    end

    if strcmp(types{t},'orf')
        cub_global_orf = cub_global;
        cub_local_orf = cub_local;
        cub_local_CAI_orf = cub_local_CAI;
        tAI_score_orf = tAI_score;
        TDR_score_orf = TDR_score;
        charge_score_orf = charge_score;
        slow_score_orf = slow_score;
    else
        cub_global_win = cub_global;
        cub_local_win = cub_local;
        cub_local_CAI_win = cub_local_CAI;
        tAI_score_win = tAI_score;
        TDR_score_win = TDR_score;
        charge_score_win = charge_score;
        slow_score_win = slow_score;
    end
end

end