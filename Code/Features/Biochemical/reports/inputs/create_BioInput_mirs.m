function [] = create_BioInput_mirs(dataset,missing_miRs)
%Creating input files for the biochemical+ model (Kd, occ features)
%dataset - A folder in miRNA/Data that contains a genes table and miRs Map Container
load(sprintf('../../../../../Data/%s/miRs.mat',dataset),'miRs');

keys = miRs.keys;
writecell(lower(replace(keys,'-','_'))',sprintf('%s/miRs.txt',dataset));

% writes general mirseqs file with all miRNAs
txt_header='mir\tguide_seq\tpass_seq\tguide_family\tpass_family';

txt=txt_header;
for i=1:length(keys)
    txt=[txt,'\n',lower(replace(keys{i},'-','_')),'\t',replace(miRs(keys{i}),'U','T'),'\t',replace(miRs(keys{i}),'U','T'),'\t',lower(replace(keys{i},'-','_')),'\t',lower(replace(keys{i},'-','_'))];  
end

fid=fopen(sprintf('%s/mirseqs.txt',dataset),'wt');
fprintf(fid,txt);
fclose(fid);

% writes mirseqs file with only miRNAs that need to be processed prior to
% Biochemical model
if ~isempty(missing_miRs)
    keys=keys(missing_miRs); % Runs feature calcs only for miRNAs that do not have files
    writecell(lower(replace(keys,'-','_'))',sprintf('%s/miRs_rel.txt',dataset));

    txt_header='mir\tguide_seq\tpass_seq\tguide_family\tpass_family';
    
    txt=txt_header;
    for i=1:length(keys)
        txt=[txt,'\n',lower(replace(keys{i},'-','_')),'\t',replace(miRs(keys{i}),'U','T'),'\t',replace(miRs(keys{i}),'U','T'),'\t',lower(replace(keys{i},'-','_')),'\t',lower(replace(keys{i},'-','_'))];  
    end

    fid=fopen(sprintf('%s/mirseqs_rel.txt',dataset),'wt');
    fprintf(fid,txt);
    fclose(fid);
end




end