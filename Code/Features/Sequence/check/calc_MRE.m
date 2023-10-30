function MREs = calc_MRE(ORF,miRNA,percents)
%Calculating number of miRNA Recognition Elements (MREs) by sequence;
%number of sites with at least p% WC/wobble base-pairing

%ORF/UTR3/miRNA - The corresponding sequences
%percents - Vector of percentages

addpath('../../Utils');

MREs = zeros(size(percents));

ORF = upper(ORF);
ORF(ORF == 'U') = 'T';
miRNA = upper(miRNA);
miRNA(miRNA == 'U') = 'T';
miRNA = miRNA(9:end);

for j = 1:length(ORF) - length(miRNA)
    seq = ORF(j:j+length(miRNA)-1);
    tail_pairing = sum(is_watson_crick(seq,miRNA,1));
    for p = 1:length(percents)
        if tail_pairing >= (percents(p) / 100) * length(miRNA)
            MREs(p) = MREs(p) + 1;
        else
            break
        end
    end
end

end