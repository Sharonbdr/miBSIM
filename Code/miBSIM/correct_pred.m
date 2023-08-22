function [predictions,features]=correct_pred(features,genes) 
% Corrects Integrative model's per site repression based on neighboring
% sites, using miBSIM.
%
% features  =  A table containing potential binding sites and their features,
%              including nominal repression S predicted by the Integrative model.
% genes     =  gene table containing gene Name, Sequence and Conservation
%
% 
% Returns   =  prediction table which include all miRNA, mRNA pairs with
%              potential sites and their predicted repression; complete feature table
%              with all binding sites and calculated features, including per site
%               repression.
%

% Variables
load('../../Data/miBSIM_vars.mat','best_c');

% Calculates DS
DS=zeros(height(features),1);
for i=1:height(features)
	DS(i)=1-sum((abs(features.S_all{i}).*exp(features.D_all{i}/(-best_c(2))))/(best_c(1)*features.site_number(i)));
end
DS(features.site_number==1)=1;

% Aggregates results and predicts corrected repression
features.miRNA = strrep(features.miRNA,'-','_');
couples=unique(features.couple_num)';

miBSIM=features.S;
predictions=[];
rep_genes={};
miRNA={};
counter=1;
for cp=couples
        gene_rows = find(features.couple_num==cp);
        not_min_TH=1;
        if length(gene_rows)>1
            %chooses validation genes that have 3'utr 7-8nt
                if any(~strcmp(features.seed_type(gene_rows),'6mer') & strcmp(features.region(gene_rows),'UTR3'))
                    curr_tab = DS(gene_rows,:).*features.S(gene_rows);
                    predictions = [predictions;[sum(features.S(gene_rows)), sum(curr_tab)]];

                    miBSIM(gene_rows)=curr_tab;
                    not_min_TH=0;
                end
             
        elseif ~strcmp(features.seed_type(gene_rows),'6mer') & strcmp(features.region(gene_rows),'UTR3')
            predictions = [predictions;[sum(features.S(gene_rows)), sum(features.S(gene_rows))]];
            not_min_TH=0;
            
        end

        if not_min_TH
            predictions = [predictions;0,0];
        end
        rep_genes{counter,1}=genes.Name{features.gene_num(gene_rows(1))};
        miRNA{counter,1}=features.miRNA{gene_rows(1)};
        counter=counter+1;

        
end

features=[features,table(DS,miBSIM,'VariableNames',{'DS','repress'})];
predictions=table(rep_genes,miRNA,predictions(:,1),predictions(:,2),'VariableNames',{'mRNA','miRNA','Integrative model','miBSIM repression'});

end
