function genes = add_cons_table(start_gene,end_gene,dataset,job_ind)
%Adding phastCons and phyloP vectors to gene table
%start/end gene - Row range in genes.mat
%dataset - A folder in miRNA/Data that contains a genes table and miRs Map Container
%job_ind - Index for the job

load(sprintf('../../Data/%s/genes.mat',dataset),'genes');

genes = genes(start_gene:end_gene,:);

for i = 1:height(genes)
    chr = genes.chr(i);
    strand = genes.strand(i);
    loc_UTR5 = genes.loc_UTR5{i};
    loc_ORF = genes.loc_ORF{i};
    loc_UTR3 = genes.loc_UTR3{i};
    [phastcons100,phylops100,phastcons20,phylops20] = add_cons(chr,strand,loc_UTR5,loc_ORF,loc_UTR3);
    genes.phastcons100(i) = {phastcons100};
    genes.phylops100(i) = {phylops100};
    genes.phastcons20(i) = {phastcons20};
    genes.phylops20(i) = {phylops20};
end

save(sprintf('../../Output/genes_%s_%d.mat',dataset,job_ind),'genes');