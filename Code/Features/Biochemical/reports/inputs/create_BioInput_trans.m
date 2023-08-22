function [] = create_BioInput_trans(dataset)
%Creating input files for the biochemical+ model (Kd, occ features)
%dataset - A folder in miRNA/Data that contains a genes table and miRs Map Container

load(sprintf('../../../../../Data/%s/genes.mat',dataset),'genes');


% transcripts
txt_header='transcript\torf\torf_length\tutr3\tutr3_length\torf_utr3';

txt=txt_header;
trans=struct;
for i = 1:height(genes)

    ORF = genes.ORF{i};
    UTR3 = genes.UTR3{i};

    ORF=replace(ORF,'U','T');
    UTR3=replace(UTR3,'U','T');

    txt=[txt,'\n',genes.Name{i},'\t',ORF,'\t',num2str(length(ORF)),'\t',UTR3,'\t',num2str(length(UTR3)),'\t',ORF,UTR3];  
    trans(i).Header=genes.Name{i};
    trans(i).Sequence=[ORF,UTR3];
end
fid=fopen(sprintf('%s/transcripts.txt',dataset),'wt');
fprintf(fid,txt);
fclose(fid);
fastawrite(sprintf('%s/orf_utr3.fa',dataset),trans);
end