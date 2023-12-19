function [min_dist, dist_start, dist_end] = calc_utr_pos(RNA_start,region_start,region_end,miRNA,seed_length)
%Calculating the relative distance from the nearest region terminus

%RNA_start = Start position of the target site in the RNA
%region_start/end = Start/End position of the region in the RNA
%miRNA = The entire miRNA strand
%seed_length = Length between 6 and 7

miR_end = RNA_start + seed_length;
miR_start = miR_end - length(miRNA) + 1;

if isinf(log10(abs(miR_start - region_start + 1)))
    dist_start = 0;
else
    dist_start = log10(abs(miR_start - region_start + 1));
end

if isinf(log10(abs(region_end - miR_end + 1)))
    dist_end = 0;
else
    dist_end = log10(abs(region_end - miR_end + 1));
end

min_dist = min(dist_start, dist_end);

end