function  sps = calc_sps(RNA,miRNA,RNA_start,seed_type,seed_length,glob_path,vienna_dir)
%Calculating SPS - The predicted seed-pairing stability, aka binding energy
%using nearest-neighbor method

%RNA = The whole RNA strand
%miRNA = The whole miRNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7
%window = Length of the upstream/downtream window before/after the RNA seed

addpath('../../Utils');

cd(vienna_dir);

%Assuming seed in miRNA starts in position 2
target_seed = RNA(RNA_start:RNA_start+seed_length-1);
if strcmp(seed_type,'8mer') || strcmp(seed_type,'7mer-A1')
    miR_seed = miRNA(1:seed_length);
else
    miR_seed = miRNA(2:seed_length+1);
end

fid = fopen('up.txt', 'wt');
fprintf(fid, '%s&%s\n@', target_seed, miR_seed);
fclose(fid);
[~,output] = system('./RNAup -o < up.txt');
output = (regexp(output,'-?\d+\.\d*','match'));
sps = str2double(output(1));

cd(glob_path);
cd '..'

end