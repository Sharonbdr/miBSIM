function  features = prepData(features)
% Adds neighboring sites features: their nominal repression S_all, the
% distance to each one [nt] D_all, and number of sites site_number.
%
% features  =  A table containing potential binding sites and their features,
%              including nominal repression S predicted by the Integrative model.
%

% Changes zeros to really small number
features.S(find(~features.S))=-min(abs(features.S(find(features.S))))/1000;

% Adds D and S to features
features = add_all_sites(features);

end



%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function features = add_all_sites(old_feat)
    
    mirs=unique(old_feat.miRNA);
    genes=unique(old_feat.gene_num)';
    features=old_feat;
    couple_counter=1;



    for g=genes
        for m=1:length(mirs)
            rows=find(old_feat.gene_num==g & strcmp(old_feat.miRNA,mirs{m}));
            N=length(rows);
            if N>1
                for i_row=1:length(rows)
                    D=old_feat.RNA_start(rows)-old_feat.RNA_start(rows(i_row));
                    S=old_feat.S(rows);
                    
                    D(i_row)=[];
                    S(i_row)=[];
               
                    features.site_number(rows(i_row))=N;
                    features.D_all{rows(i_row)}=abs(D);
                    features.S_all{rows(i_row)}=S;
                    features.couple_num(rows(i_row))=couple_counter;
                end
            elseif N==1
                features.site_number(rows)=N;
                features.D_all{rows}=inf;
                features.S_all{rows}=0;
                features.couple_num(rows)=couple_counter;
            end

           couple_counter=couple_counter+1;

        end
    end
end
