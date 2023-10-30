function [features] = add_bio_features(genes,features,cnn_results_path)
% Retrieves biochemical features calculated by the Biochemical+ model, and
% adds them to the base feature set.
%
% genes            =  gene table containing gene Name, Sequence and Conservation
%                     variables
% features         =  A table containing potential binding sites and their
%                     features, excluding the biochemical features
% cnn_results_path =  str of the path to CNN prediction file ending with miRNA.txt
%
% Returns          =  features table with added biochemical features Kd occ and occ_init 
%

% adds Kd, occ, occ_init features
mirs=replace(unique(features.miRNA)','-','_');
gene_num = containers.Map(1:height(genes),genes.Name);
gene_name = containers.Map(genes.Name,1:height(genes));


for m=mirs
    rep_cnn=readtable([cnn_results_path,sprintf('%s.txt',m{1})]);
    rows=find(strcmp(rep_cnn{:,8},'no site'));
    rep_cnn(rows,:)=[];

    n=width(rep_cnn);
    if n<24
        error('Biochemical model returned NaN prediction, please rerun the biochemical features')
    end

    val=rep_cnn{1,7}; % matlab sometimes reads the column names correctly and sometimes not, this adjusts for both
    if ~iscell(val) % read correctly
        rows=find(strcmp(rep_cnn{:,8},'no site'));
    else            % did not read correctly
        rows=find(strcmp(rep_cnn{:,7},'no'));
    end
    rep_cnn(rows,:)=[];
    
    % matches the sites by location on mRNA strand
    rows=find(strcmp(features.miRNA,m{1}));
    temp=features(rows,:);
    for i=1:height(temp)
        curr_g=temp.gene_num(i);
        cnn_rows=find(strcmp(gene_num(curr_g),rep_cnn{:,2}));
        global_ind=rows(i);
        if ~isempty(cnn_rows)
            temp_cnn=rep_cnn(cnn_rows,:);
    
            [~,site_row]=min(abs(temp_cnn{:,5}-temp.RNA_start(i)));

            features.kds_val(global_ind)=exp(temp_cnn{site_row,4});
            features.occ(global_ind)=temp_cnn{site_row,end-1};
            features.occ_init(global_ind)=temp_cnn{site_row,end};
            features.logKd(global_ind)=temp_cnn{site_row,4};
            features.log_occ(global_ind)=log(abs(temp_cnn{site_row,end-1}));
            features.log_occ_init(global_ind)=log(abs(temp_cnn{site_row,end}));
        else
            features.kds_val(global_ind)=NaN;
            features.occ(global_ind)=NaN;
            features.occ_init(global_ind)=NaN;
            features.logKd(global_ind)=NaN;
            features.log_occ(global_ind)=NaN;
            features.log_occ_init(global_ind)=NaN;
        end
    end
end

end