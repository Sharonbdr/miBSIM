function score = align3p(RNA_3p_seq,miR_3p_seq,miR_offset,RNA_offset,overhang)
%Helper function for calc_3p_pairing - calculating pairing score
score = 0;
tempscore = 0;
prevmatch = 0;
i = 0;
offset = max(miR_offset,RNA_offset);
addpath('../../Utils')

while (i<length(miR_3p_seq)-miR_offset) && (i<length(RNA_3p_seq)-RNA_offset)
    if is_watson_crick(RNA_3p_seq(i+1+RNA_offset),miR_3p_seq(i+1+miR_offset))
       if (i+miR_offset-overhang>=4)  && (i+miR_offset-overhang<=7)
           if prevmatch == 0
               tempscore = 0;
           end
           tempscore = tempscore + 1;
       else
           if prevmatch == 0
               tempscore = 0;
           end
           tempscore = tempscore + 0.5;
       end
       prevmatch = prevmatch + 1;
    elseif prevmatch >= 2
        if tempscore > score
            score = tempscore;
        end
        tempscore = 0;
        prevmatch = 0;
    else
        tempscore = 0;
        prevmatch = 0;
    end
    i = i + 1;
end

if prevmatch >=2 && tempscore > score
    score = tempscore;
end

score = score - max(0,((offset-2)/2));