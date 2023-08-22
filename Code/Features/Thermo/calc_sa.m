function  sa = calc_sa(RNA,RNA_start,seed_type,region_end,region,glob_path,vienna_dir)
%Calculating SA - The structural accessibility of the target site

%RNA = The whole RNA strand
%miRNA = The whole miRNA strand
%region_end = End position of the region in the RNA
%region = 5'UTR/ORF/3'UTR

addpath('../../Utils');

cd(vienna_dir);

if any(strcmp(seed_type,{'8mer','7mer-m8','offset-6mer'}))
    seven_pos = RNA_start + 1;
elseif strcmp(seed_type,'offset-7mer')
    seven_pos = RNA_start + 2;
elseif strcmp(seed_type,'6mer-A1')
    seven_pos = RNA_start - 1;
else
    seven_pos = RNA_start;
end

if not(seven_pos+6 > region_end & strcmp(region,'UTR3'))
    fid = fopen('pl.txt', 'wt');
    fprintf(fid, '%s\n@', RNA);
    fclose(fid);
    system('./RNAplfold -L 40 -W 80 -u 14 < pl.txt');
    fid = fopen('plfold_lunp');
    for j = 1:seven_pos+6+2
        line = fgetl(fid);
    end
    fclose(fid);
    prob = str2double(regexp(line,'0\.\d*','match'));
    if isempty(prob)
        sa = 0;
    else
        sa = log10(prob(end));
    end
else
    sa = NaN;
end

cd(glob_path);
cd '..';

end