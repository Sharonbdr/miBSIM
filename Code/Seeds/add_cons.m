function [phastcons100,phylops100,phastcons20,phylops20] = add_cons(chr,strand,loc_UTR5,loc_ORF,loc_UTR3)
%Extracting phastCons and phyloP vectors from hg38 chromosomal files
%chr - Chromomsome number (X = 23, Y = 24)
%strand - 1 if (+), -1 if (-)
%loc_UTR5/ORF/UTR3 - Genomic coordinates (start and end) of regions

locs_array = cell(1,3);
locs_array{1} = loc_UTR5;
locs_array{2} = loc_ORF;
locs_array{3} = loc_UTR3;

gene_coordinates = [];

for r = 1:3
    locs = locs_array{r};
    coordinates = [];
    for j = 1:size(locs,1)
        if strand == 1
            start_coor = locs(j,1);
            end_coor = locs(j,2);
            coordinates = [coordinates start_coor:end_coor]; %#ok<*AGROW>
        %Coordinates are reversed
        else
            start_coor = locs(j,2);
            end_coor = locs(j,1);
            coordinates = [coordinates start_coor:-1:end_coor];
        end
    end         
    gene_coordinates = [gene_coordinates coordinates];
end

load(sprintf('../../Data/Conservation/phastCons/phastcons%d.mat',chr));
phastcons100 = phastcons(gene_coordinates);
clear phastcons
load(sprintf('../../Data/Conservation/phylop/phylops%d.mat',chr));
phylops100 = phylops(gene_coordinates);
clear phylops
load(sprintf('../../Data/Conservation/phastCons20/phastcons%d.mat',chr));
phastcons20 = phastcons(gene_coordinates);
clear phastcons
load(sprintf('../../Data/Conservation/phylop20/phylops%d.mat',chr));
phylops20 = phylops(gene_coordinates);
clear phylops

end