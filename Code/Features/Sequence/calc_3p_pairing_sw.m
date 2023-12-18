function pairing3p_sw = calc_3p_pairing_sw(RNA,RNA_start,miRNA,seed_type)
%Calculating the 3' complementary pairing score

%RNA = The whole RNA strand
%RNA_start = Start position of the seed in the RNA
%miRNA = The whole miRNA strand
%seed_length = Length between 6 and 7

addpath('../../Utils');

sw = 7; %Length of sliding window

if strcmp(seed_type,'8mer') || strcmp(seed_type,'7mer-m8')
    miR_3p_len = length(miRNA) - 8;
    miR_3p_1 = 9; %The first position after the seed
else
    miR_3p_len = length(miRNA) - 7;
    miR_3p_1 = 8;
end

RNA_5p_seq = RNA(max(1,RNA_start - miR_3p_len):RNA_start - 1);

pairing3p_sw = zeros(1,length(RNA_5p_seq)-sw);
for i = 0:length(RNA_5p_seq)-sw-1
    pairing3p_sw(i+1) = sum(is_watson_crick(miRNA(miR_3p_1+i:miR_3p_1+i+sw-1),RNA_5p_seq(i+1:i+sw)));
end

pairing3p_sw = mean(pairing3p_sw);