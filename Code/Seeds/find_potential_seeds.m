function seed_table = find_potential_seeds(RNA,miR,UTR5_end,ORF_end,is_canonical,keep_all)
% Find potential seeds in RNA according to miRNA strand
%
% RNA             =  The whole RNA strand
% miR             =  The whole miRNA strand
% UTR5/ORF_end    =  End position of the region in the RNA
% is_canonical    =  'canonical' for only canonical sites, 'all' for all sites
% keep_all        =  whether or not we should keep all possible sites, or include
%                    only 3'UTR sites and sites that reside in the last third of the orf
%

if ~strcmp(is_canonical,'canonical') && ~strcmp(is_canonical,'all')
    error('Please enter ''canonical'' or ''all''')
end

RNA = upper(RNA);
RNA(RNA == 'U') = 'T';
miR = upper(miR);
miR(miR == 'U') = 'T';

seed_table = table(0);

sites_8mer = regexp(RNA,[seqrcomplement(miR(8)),'(?=',seqrcomplement(miR(2:7)),'A)'])';
if ~isempty(sites_8mer)
    temp_table = table(sites_8mer,cellstr(repmat('8mer',length(sites_8mer),1)),repmat(8,length(sites_8mer),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    seed_table = temp_table;
end

sites_7mer_m8 = regexp(RNA,[seqrcomplement(miR(8)),'(?=',seqrcomplement(miR(2:7)),'[^A])'])';
if ~isempty(sites_7mer_m8)
    temp_table = table(sites_7mer_m8,cellstr(repmat('7mer-m8',length(sites_7mer_m8),1)),repmat(7,length(sites_7mer_m8),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

sites_7mer_A1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miR(8))),seqrcomplement(miR(2:7)),'A)'])';
if ~isempty(sites_7mer_A1)
    temp_table = table(sites_7mer_A1+1,cellstr(repmat('7mer-A1',length(sites_7mer_A1),1)),repmat(7,length(sites_7mer_A1),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

sites_6mer = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miR(8))),seqrcomplement(miR(2:7)),'[^A])'])';
if ~isempty(sites_6mer)
    temp_table = table(sites_6mer+1,cellstr(repmat('6mer',length(sites_6mer),1)),repmat(6,length(sites_6mer),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

if strcmp(is_canonical,'all')
    sites_6mer_A1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miR(7))),seqrcomplement(miR(2:6)),'A)'])';
    if ~isempty(sites_6mer_A1)
        temp_table = table(sites_6mer_A1+1,cellstr(repmat('6mer-A1',length(sites_6mer_A1),1)),repmat(5,length(sites_6mer_A1),1),...
            'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_offset_7mer = regexp(RNA,[seqrcomplement(miR(9)),'(?=',...
        seqrcomplement(miR(3:8)),sprintf('[^%s]',seqrcomplement(miR(2))),')'])';
    if ~isempty(sites_offset_7mer)
        temp_table = table(sites_offset_7mer,cellstr(repmat('offset-7mer',length(sites_offset_7mer),1)),...
            repmat(7,length(sites_offset_7mer),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table; %#ok<*NASGU>
        else
            seed_table = [seed_table; temp_table];
        end
    end


    sites_offset_6mer = regexp(RNA,[sprintf('[^%s](?=',...
        seqrcomplement(miR(9))),seqrcomplement(miR(3:8)),sprintf('[^%s]',seqrcomplement(miR(2))),')'])';
    if ~isempty(sites_offset_6mer)
        temp_table = table(sites_offset_6mer+1,cellstr(repmat('offset-6mer',length(sites_offset_6mer),1)),...
            repmat(6,length(sites_offset_6mer),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_centered = regexp(RNA,[seqrcomplement(miR(13)),'(?=',seqrcomplement(miR(3:12)),')|',seqrcomplement(miR(14)),'(?=',seqrcomplement(miR(4:13)),')|',seqrcomplement(miR(15)),'(?=',seqrcomplement(miR(5:14)),')'])';
    if ~isempty(sites_centered)
        temp_table = table(sites_centered+1,cellstr(repmat('centered',length(sites_centered),1)),...
            repmat(5,length(sites_centered),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end
    
    sites_CDNST1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miR(7))),seqrcomplement(miR(2:6)),'[^A])'])';
    if ~isempty(sites_CDNST1)
        temp_table = table(sites_CDNST1+1,cellstr(repmat('CDNST1',length(sites_CDNST1),1)),...
            repmat(5,length(sites_CDNST1),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST2 = regexp(RNA,[seqrcomplement(miR(7)),'(?=',...
        seqrcomplement(miR(6)),sprintf('[^%s]',seqrcomplement(miR(5))),seqrcomplement(miR(2:4)),'A)'])';
    if ~isempty(sites_CDNST2)
        temp_table = table(sites_CDNST2,cellstr(repmat('CDNST2',length(sites_CDNST2),1)),...
            repmat(6,length(sites_CDNST2),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST3 = regexp(RNA,[seqrcomplement(miR(9)),'(?=',seqrcomplement(miR(7:8)),sprintf('[^%s]',seqrcomplement(miR(6))),seqrcomplement(miR(5)),...
        sprintf('[^%s]',seqrcomplement(miR(4))),seqrcomplement(miR(3)),sprintf('[^%s]',seqrcomplement(miR(2))),')'])';
    if ~isempty(sites_CDNST3)
        temp_table = table(sites_CDNST3,cellstr(repmat('CDNST3',length(sites_CDNST3),1)),...
            repmat(7,length(sites_CDNST3),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST4 = regexp(RNA,[seqrcomplement(miR(8)),'(?=',sprintf('[^%s]',seqrcomplement(miR(7))),...
        sprintf('[^%s]',seqrcomplement(miR(6))),sprintf('[^%s]',seqrcomplement(miR(5))),seqrcomplement(miR(2:4)),')'])';
    if ~isempty(sites_CDNST4)
        temp_table = table(sites_CDNST4,cellstr(repmat('CDNST4',length(sites_CDNST4),1)),...
            repmat(3,length(sites_CDNST4),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end
end


if width(seed_table) == 1
    seed_table(1,:) = [];
else
    seed_table{seed_table.RNA_start <= UTR5_end,'region'} = {'UTR5'};
    seed_table{seed_table.RNA_start > UTR5_end & seed_table.RNA_start <= ORF_end ,'region'} = {'ORF'};
    seed_table{seed_table.RNA_start > ORF_end,'region'} = {'UTR3'};

    if ~keep_all
        seed_table(seed_table.RNA_start < ORF_end*2/3,:) =[];
    end
end



end