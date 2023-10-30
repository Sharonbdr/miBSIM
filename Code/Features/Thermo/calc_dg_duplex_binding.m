function  [dg_duplex, dg_duplex_seed, dg_binding, dg_binding_seed] = ...
    calc_dg_duplex_binding(RNA,miRNA,miR_start,RNA_pos1,RNA_start,seed_type,seed_length,glob_path,vienna_dir)
%Calculating dg_duplex - The MFE of the miRNA-RNA complex; and dg_binding -
%The binding energy of the RNA and miRNA

%RNA = The whole RNA strand
%miRNA = The whole miRNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7

addpath('../../Utils');

cd(vienna_dir);

if strcmp(seed_type,'8mer') || strcmp(seed_type,'7mer-A1')
    pairing = is_watson_crick(RNA(RNA_start:RNA_start+seed_length-2), miRNA(miR_start:miR_start+seed_length-2));
else
    pairing = is_watson_crick(RNA(RNA_start:RNA_start+seed_length-1), miRNA(miR_start:miR_start+seed_length-1));
end

seed_constraints = '';
site_constraints = '';

%Constraints: RNA & miRNA match on seed (If Watson-Crick), no constraints
%on other nts
for i = 1:length(pairing)
    if ~pairing(i)
        site_constraints = [site_constraints, '.']; %#ok<*AGROW>
    else
        site_constraints = [site_constraints, '('];
    end
end

for i = length(pairing):-1:1
    if ~pairing(i)
        seed_constraints = [seed_constraints, '.'];
    else
        seed_constraints = [seed_constraints, ')'];
    end
end

if strcmp(seed_type,'8mer') || strcmp(seed_type,'7mer-m8')
    len_no_constraints = length(miRNA) - 8;
else
    len_no_constraints = length(miRNA) - 7;
end

%Input parameters for RNAcofold

start_theoretic = RNA_pos1 - (length(miRNA)-1);
end_theoretic = RNA_pos1;

%If the window overhangs from the 3' end, add polyA tail
%If the miRNA overhangs from the 5'UTR, take only the overlapping part
if end_theoretic > length(RNA)
    disp('Window exceeds UTR length, adding polyA tail');
    target_site = RNA(max(1,start_theoretic):end);
    target_site = [target_site, repmat('A',1,end_theoretic-length(RNA))];
else
    target_site = RNA(max(1,start_theoretic):end_theoretic);
end

constraints = [repmat('.',1,min(len_no_constraints,RNA_start - 1)), site_constraints, ...
 '.&.',seed_constraints,repmat('.',1,len_no_constraints)];

fid = fopen('cofold.txt', 'wt');
fprintf(fid, '%s&%s\n%s\n@', target_site, miRNA, constraints);
fclose(fid);

%target_site and miRNA
[~,output] = system('./RNAcofold --constraint --partfunc < cofold.txt'); disp(output)
output = (regexp(output,'-?\d+\.\d*','match'));
dg_duplex = str2double(output(1));
dg_binding = str2double(output(end));

%target_seed and miR_seed
target_seed = RNA(RNA_start:RNA_start+seed_length-1);
miR_seed = miRNA(miR_start:miR_start+seed_length-1);

fid = fopen('cofold.txt', 'wt');
fprintf(fid, '%s&%s\n@', target_seed, miR_seed);
fclose(fid);
[~,output] = system('./RNAcofold --partfunc < cofold.txt');
output = (regexp(output,'-?\d+\.\d*','match'));
dg_duplex_seed = str2double(output(1));
dg_binding_seed = str2double(output(end));

cd(glob_path);
cd('..')