function MRE = calc_MRE(ORF,miRNA,max_energy,mismatch,glob_path,vienna_dir)
%Counting number of MREs - sites with sufficiently low binding energy

%ORF = The whole ORF
%miRNA = The whole miRNA strand
%max_energy = Vector of maximal accepted energies for MRE
cd(vienna_dir);

%Input parameters for RNAduplex

fid = fopen('duplex.txt', 'wt');
fprintf(fid, '%s\n%s\n@', ORF, miRNA);
fclose(fid);

%target_site and miRNA
% vienna_dir
[~,output] = system('./RNAduplex -e 100 < duplex.txt');
if ~isempty(output)
output = splitlines(output);
output = output(contains(output,'&'));
output = extractAfter(output,'&');
output = table(output,'VariableNames',{'string'});

for i = 1:height(output)
    str = output.string{i};
    %Energy
    val = regexp(str,'\( ?\-?\d+\.\d+\)','match');
    output.energy(i) = str2double(extractBetween(val,'(',')','Boundaries','exclusive'));
    %Start position
    val = regexp(str,'\d+,','match');
    start_pos = str2num(val{2});
    %End position
    val = regexp(str,',\d+','match');
    end_pos = str2num(val{2});
    %Pairing
    output.pairing{i} = extractBefore(str,' ');
    
    if start_pos <= 8
        pairing = repmat('.',1,length(miRNA));
        pairing(start_pos:end_pos) = output.pairing{i};
        if sum(pairing(2:8) == '.') < mismatch
            output.energy(i) = 999;
        end
    end
end

MRE = [];
diff = abs(max_energy(2) - max_energy(1));
for m = max_energy
    MRE = [MRE,sum(output.energy <= m & output.energy > m-diff)]; %#ok<AGROW>
end

cd(glob_path);
cd '..'
else
    MRE=zeros(size(max_energy));
    cd(glob_path);
    cd '..'
end
end