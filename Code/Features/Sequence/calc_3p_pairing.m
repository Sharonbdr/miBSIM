function pairing3p = calc_3p_pairing(RNA,RNA_start,miRNA,seed_length,miR_start)
%Calculating the 3' complementary pairing score: +1 pt for match to nt
%13-16 of miRNA; +0.5 pt for extending the match; -0.5 pt for every offset
%nt beyond 2

%RNA = The whole RNA strand
%RNA_start = Start position of the seed in the RNA
%miRNA = The whole miRNA strand
%seed_length = Length of seed
%miR_start = Beginning position of seed in miRNA

addpath('../../Utils');

pairing3p = 0;

seed_end = miR_start + seed_length - 1;
min_offset = seed_end - 13 + 1;
max_offset = length(miRNA) - 16;
RNA_nt1 = RNA_start+seed_length-1+miR_start-1; %The position across from nt1 of miRNA

%Adding Poly A tail in case the miRNA protrudes from 3' end
added_A = RNA_nt1 - length(RNA);
RNA = [RNA, repmat('A',1,added_A)];

for offset = min_offset:max_offset
    if RNA_nt1-15-offset >= 1
        score = sum(is_watson_crick(miRNA(13:16),RNA(RNA_nt1-15-offset:RNA_nt1-12-offset))) - max(0,((abs(offset)-2)/2));
        i = 1;
        while 13-i > seed_end && is_watson_crick(miRNA(13-(i-1)),RNA(RNA_nt1-12-offset+(i-1)))
            if is_watson_crick(miRNA(13-i),RNA(RNA_nt1-12-offset+i))
                score = score + 0.5;
            end
            i = i + 1;
        end

        i = 1;
        while 16+i <= length(miRNA) && RNA_nt1-15-offset-i >= 1 && is_watson_crick(miRNA(16+(i-1)),RNA(RNA_nt1-15-offset-(i-1)))
            if is_watson_crick(miRNA(16+i),RNA(RNA_nt1-15-offset-i))
                score = score + 0.5;
            end
            i = i + 1;
        end

        if score > pairing3p
            pairing3p = score;
        end
    end
end

end