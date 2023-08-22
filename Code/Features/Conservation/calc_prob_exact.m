function prob_exact = calc_prob_exact(RNA, RNA_start, seed_length)
%Calculating the exact probability of the target site

%RNA = The whole RNA strand
%RNA_start = Start position of the seed
%seed_length = Length between 6 and 7

RNA = lower(RNA);
RNA(RNA == 't') = 'u';
target_site = RNA(RNA_start:RNA_start + seed_length - 1);
nobs = count(RNA,target_site);
alphabet = 'acgu';
[~,trans_matrix] = calc_transition_dict(RNA);

fname = 'transitions.markov';
fid = fopen(fname,'wt');
trans_str = '';

for i = 1:length(alphabet)
    for j = 1:length(alphabet)
        trans_str = [trans_str, ' ', num2str(trans_matrix(i,j))];
    end
    trans_str = [trans_str, '\n'];
end

fprintf(fid,trans_str);
fclose(fid);

command = sprintf('./spatt --pattern %s --nobs %d --seqlen %d --alphabet %s -m 1 --over --Markov %s',...
    target_site, nobs, length(RNA), alphabet,fname);
[~,output] = system(command);
index = regexp(output,'=');
index = index(end);
prob_exact = str2double(extractAfter(output,index));

delete(fname);

end