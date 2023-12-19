function permuted_ORF = permute_orf(ORF)
%Permute ORF, while keeping nucleotide and amino acid composition

ORF = upper(ORF);
ORF(ORF == 'U') = 'T';

if mod(length(ORF),3) ~= 0
    error('ORF length should be divisible by 3');
end

load('../../Data/codon_table.mat')
aas = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

orf_codons = regexp(ORF,'.{3}','match');

for i = 1:length(aas)
    curr_aa = aas{i};
    codon_inds = [];
    if ~any(strcmp(curr_aa,{'M','W'})) %Since M and W have only 1 codon
        curr_codons = codon_table(curr_aa);
        for j = 1:length(curr_codons)
            codon_inds = [codon_inds,find(strcmp(orf_codons,curr_codons{j}))]; %#ok<AGROW>
        end
        permutation = randperm(length(codon_inds));
        permuted_codon_inds = codon_inds(permutation);
        orf_codons(codon_inds) = orf_codons(permuted_codon_inds);
    end
end

permuted_ORF = strjoin(orf_codons,'');

end