function Input = create_BioInput(dataset)
%Creating input files for the biochemical+ model (Kd, occ features)
%dataset - A folder in miRNA/Data that contains a genes table and miRs Map Container


load(sprintf('../../Data/%s/genes.mat',dataset),'genes');
load(sprintf('../../Data/%s/miRs.mat',dataset),'miRs');

Input = table(0);
keys = miRs.keys;

% miRs
txt_header='mir\tguide_seq\tpass_seq\tguide_family\tpass_family';

txt=txt_header;
for i=1:length(keys)
    txt=[txt,'\n',replace(keys{i},'-','_'),'\t',replace(miRs(keys{i}),'U','T'),'\t',replace(miRs(keys{i}),'U','T'),'\t',replace(keys{i},'-','_'),'\t',replace(keys{i},'-','_')];  
end

fid=fopen('../Biochemical/reports/inputs/mirseqs.txt','wt');
fprintf(fid,txt);
fclose(fid);

% transcripts
txt_header='transcript\torf\torf_length\tutr3\tutr3_length\torf_utr3';

txt=txt_header;
trans=Table;
for i = 1:length(genes)

    ORF = genes.ORF{i};
    UTR3 = genes.UTR3{i};

    ORF=replace(ORF,'U','T');
    UTR3=replace(UTR3,'U','T');

    txt=[txt,'\n',genes.Name{i},'\t',ORF,'\t',length(ORF),'\t',UTR3,'\t',length(UTR3),'\t',ORF,UTR3];  
    trans.Header{i}=genes.Name{i};
    trans.Sequence{i}=[ORF,UTR3];
end
fid=fopen('../Biochemical/reports/inputs/transcripts.txt','wt');
fprintf(fid,txt);
fclose(fid);
fastawrite('../Biochemical/reports/inputs/orf_utr3.fa',trans);
end