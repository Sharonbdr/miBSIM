function emp_p = calc_emp_p(RNA,RNA_start,seed_length,ORF_start,ORF_end)
%Calculating the empirical p-value of the site, randomizing the ORF
%(while keeping the primary structure)

%RNA = The whole RNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7
%ORF_start/end - Start position of the ORF in the gene

cd('../../Utils')

perm_num = 10000;
RNA(RNA == 'U') = 'T';

ORF_seed_start = RNA_start - ORF_start + 1;
ORF_seed_end = min(ORF_seed_start + seed_length-1,ORF_end - ORF_start + 1);
ORF = RNA(ORF_start:ORF_end);
seed = ORF(ORF_seed_start:ORF_seed_end);

    function seed_exists = find_seed(random_ORF)
        seed_exists = strcmp(seed,random_ORF(ORF_seed_start:ORF_seed_end));
    end

orf_array = cell(perm_num,1);
    orf_array(:) = {ORF};
    orf_array = cellfun(@permute_orf,orf_array,'UniformOutput',0);
    emp_p = mean(cellfun(@find_seed,orf_array));

cd('../Features/Conservation')

end