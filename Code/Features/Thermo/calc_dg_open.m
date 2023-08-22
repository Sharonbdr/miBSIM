function dg_open = calc_dg_open(RNA,RNA_pos1,miRNA_length,glob_path,vienna_dir)
%Calculating dg_open - The energy needed to maintain the RNA open in order to
%accomodate for RISC

%RNA = The whole RNA strand
%RNA_start = Start position of the seed
%seed_length = Length between 6 and 7

%Default parameters - Length of RISC up/downstream of seed (that needs to stay open)
RISC_win = floor((70 - miRNA_length) / 2);

addpath('../../Utils');

cd(vienna_dir);

start_theoretic = RNA_pos1 - miRNA_length - RISC_win;
start = max(1,start_theoretic);
%If RISC_win exceeds the 3' end, add polyA tail
end_theoretic = RNA_pos1 + RISC_win;
if end_theoretic > length(RNA)
    curr_RNA = RNA(start:end);
    curr_RNA = [curr_RNA, repmat('A',1,end_theoretic-length(RNA))];
else
    curr_RNA = RNA(start:end_theoretic);
end

%Without constraints
fid = fopen('fold.txt', 'wt');
fprintf(fid, '%s\n@', curr_RNA);
fclose(fid);
[~, output] = system('./RNAfold --partfunc < fold.txt');
dg_unconstrained = str2double(regexp(output,'-?\d+\.\d*','match'));
dg_unconstrained = dg_unconstrained(2);

%With constraints
fid = fopen('fold.txt', 'wt');
fprintf(fid, '%s\n%s\n@', curr_RNA, repmat('x',1,length(curr_RNA)));
fclose(fid);
[~, output] = system('./RNAfold -C --partfunc < fold.txt');
dg_constrained = str2double(regexp(output,'-?\d+\.\d*','match'));
dg_constrained = dg_constrained(2);

dg_open = dg_constrained - dg_unconstrained;

cd(glob_path);
cd('..')

end