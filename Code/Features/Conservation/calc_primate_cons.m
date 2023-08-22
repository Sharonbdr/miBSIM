function [changeN1,changeN2,changeN3,changeA1,changeA2,changeA3] =...
    calc_primate_cons(ENST,mut_ind,region,ORF_start,RNA_start,seed_length,dataset) %#ok<STOUT>
%Calculating the average number of differences between human and subset_i,
%on the nt and AA level
%subset_1 = {Chimp}; subset_2 = {Chimp, Gorilla, Orangutan}; subset_3 = All
%primates

%Parameters:
%name - Gene Name
%region - ORF/UTR3
%region_start - Position of the 1st nt in region
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7

load(sprintf('../../../Data/Conservation/MSA/%s/%s_%d.mat',dataset,ENST,mut_ind),'orf_msa','utr_msa');

msa = [orf_msa,utr_msa];

RNA_start_local = RNA_start - ORF_start + 1; %RNA_start in terms of region coordinates (instead of the whole gene)
RNA_end_local = RNA_start_local + seed_length - 1;

ind_vec = find(msa(1,:) ~= '-');
site_coords = ind_vec(RNA_start_local:RNA_end_local);

msa = msa(:,site_coords);

subsets = {1:2,1:4,1:12}; %subset_1-3

%Comparing nt by nt; indel = 1, transversion = 0.5, transition = 0.25

for s = 1:length(subsets)
    curr_msa = msa(subsets{s},:);
    changes = zeros(1,length(subsets{s})-1);
    for i = 2:size(curr_msa,1)
        for j = 1:seed_length
            if curr_msa(i,j) == '-'
                changes(s) = changes(s) + 1;
            elseif curr_msa(1,j) ~= curr_msa(i,j)
                if any(strcmp(curr_msa(1,j),{'A','G'})) && any(strcmp(curr_msa(i,j),{'A','G'}))
                    changes(s) = changes(s) + 0.25;
                elseif any(strcmp(curr_msa(1,j),{'C','T'})) && any(strcmp(curr_msa(i,j),{'C','T'}))
                    changes(s) = changes(s) + 0.25;
                else
                    changes(s) = changes(s) + 0.5;
                end
            end
        end
    end
    changes = mean(changes); %#ok<NASGU>
    eval(sprintf('changeN%d = changes;',s))
end

%Comparing AA by AA; penalty according to PAM30

if strcmp(region,'ORF')
    ind_vec = find(orf_msa(1,:) ~= '-');
    %Taking the whole codon
    switch mod(RNA_start_local,3)
        case 0
            msa_start = RNA_start_local - 2;
        case 1
            msa_start = RNA_start_local;
        case 2
            msa_start = RNA_start_local - 1;
    end
    msa_start = max(msa_start,1);
    
    switch mod(RNA_end_local,3)
        case 0
            msa_end = RNA_end_local;
        case 1
            msa_end = RNA_end_local + 2;
        case 2
            msa_end = RNA_end_local + 1;
    end
    msa_end = min(msa_end,length(ind_vec));
    
    site_coords = ind_vec(msa_start:msa_end);
    msa = orf_msa(:,site_coords);

    load('../../../Data/PAM.mat','PAM');
    for s = 1:length(subsets)
        curr_msa = msa(subsets{s},:);
        changes = zeros(1,length(subsets{s})-1);
        for i = 2:size(curr_msa,1)
            for j = 1:3:length(site_coords)-2
                hs_codon = upper(curr_msa(1,j:j+2));
                pri_codon = upper(curr_msa(i,j:j+2));
                codons = {hs_codon,pri_codon};
                for c = 1:length(codons)
                    if contains(codons{c},'-')
                        codons{c} = 'indel';
                    elseif nt2aa(codons{c},'AlternativeStartCodons',false,'ACGTOnly',false) == '*'
                        codons{c} = 'STOP';
                    else
                        codons{c} = nt2aa(codons{c},'AlternativeStartCodons',false,'ACGTOnly',false);
                    end
                end
                changes(i-1) = changes(i-1) + PAM{codons{1},codons{2}};
            end
        end
        changes = mean(changes); %#ok<NASGU>
        eval(sprintf('changeA%d = changes;',s))
    end

else
    changeA1 = 999;
    changeA2 = 999;
    changeA3 = 999;
end

end