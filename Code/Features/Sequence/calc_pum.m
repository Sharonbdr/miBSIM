function [ORF_pum,UTR3_pum,nei_pum] = calc_pum(RNA,RNA_start,seed_length,ORF_start,ORF_end)
%finding PUM sites in ORF and 3'UTR

motifs{1} = 'UGUA[ACGU]AU[AU]';
motifs{2} = '[CGU][CU]C[CG][AU]C[CG]C';

ORF_pum = zeros(1,length(motifs));
UTR3_pum = zeros(1,length(motifs));
nei_pum = zeros(1,length(motifs));

for i = 1:length(motifs)
    motif_inds = regexp(RNA,[motifs{i}(1),'(?=',motifs{i}(2:end),')']);
    if ~isempty(motif_inds)
        motif_inds = motif_inds(motif_inds <= RNA_start - 8 | motif_inds > RNA_start + seed_length-1);

        %1 pt for the first site, 0.5 pts for every additional site

        ORF_num = sum(motif_inds >= ORF_start & motif_inds <= ORF_end - 8);
        if ORF_num > 0
            ORF_pum(i) = 0.5*(ORF_num - 1) + 1; %#ok<*AGROW>
        else
            ORF_pum(i) = 0;
        end

        UTR3_num = sum(motif_inds >= ORF_end+1 & motif_inds <= length(RNA) - 8);
        if UTR3_num > 0
            UTR3_pum(i) = 0.5*(UTR3_num - 1) + 1; %#ok<*AGROW>
        else
            UTR3_pum(i) = 0;
        end
        
        nei_pum(i) = sum(motif_inds >= RNA_start - 400 & motif_inds + 7 <= RNA_start + seed_length-1 + 400);
    end
end
    
end