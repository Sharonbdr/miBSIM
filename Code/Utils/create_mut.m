function mut_gene = create_mut(gene,mut_ind,mut_type,mut_pos,new_nt)
%gene is a 1xn table with the usual columns: chr, UTR5/ORF/UTR3, strand, phastcons/phylops, mut_ind == 0
%mut_ind = The mutation's index (original gene = 0); an integer
%mut_type = insert/del/snp
%mut_pos = The 1xn position vector of the mutation (in case of insert - the
%positions in which the new_nt enter)
%new_nt = The 1xn char array entered in mut_pos. In case of deletion mutation, no need to input new_nt

%Making sure mut_type and new_nt are the right inputs
if strcmp(mut_type,'insert') || strcmp(mut_type,'snp')
    if nargin < 6
        error('Please input the new nt')
    end
elseif ~strcmp(mut_type,'del')
    error('Mutation type should be insert, del or snp');
end

%Checking mut_pos is a positive integer
if any(round(mut_pos) ~= mut_pos) || any(mut_pos <= 0)
    error('Position should be a positive integer')
end

%Checking the length of mut_pos and new_nt is the same
if nargin == 6 && length(mut_pos) ~= length(new_nt)
    error('The length of mut_pos and new_nt shuld be identical')
end

ORF_start = length(gene.UTR5{1}) + 1;
UTR3_start = length(gene.UTR5{1}) + length(gene.ORF{1}) + 1;

%Checking for multiple regions
all_mut_regions = cell(size(mut_pos));
mut_pos_reg = zeros(size(mut_pos)); %mut_pos relative to region start
for i = 1:length(mut_pos)
    if mut_pos(i) < ORF_start
        all_mut_regions{i} = 'UTR5';
        mut_pos_reg(i) = mut_pos(i);
    elseif mut_pos(i) < UTR3_start
        all_mut_regions{i} = 'ORF';
        mut_pos_reg(i) = mut_pos(i) - length(gene.UTR5{1});
    else
        all_mut_regions{i} = 'UTR3';
        mut_pos_reg(i) = mut_pos(i) - length(gene.UTR5{1}) - length(gene.ORF{1});
    end
end
mut_regions = unique(all_mut_regions);

%Editing sequence region-by-region
for r = 1:length(mut_regions)
    inds = find(strcmp(all_mut_regions,mut_regions{r}));
    curr_mut_pos_reg = mut_pos_reg(inds);
    if ~strcmp(mut_type,'del')
        curr_new_nt = new_nt(inds);
    end

    region_seq = gene{1,mut_regions{r}}{1};
    if strcmp(mut_type,'insert')
        region_seq = [region_seq(1:min(curr_mut_pos_reg)-1), curr_new_nt, region_seq(:,min(curr_mut_pos_reg):end)];
    elseif strcmp(mut_type,'del')
        region_seq(curr_mut_pos_reg) = [];
    elseif strcmp(mut_type,'snp')
        region_seq(curr_mut_pos_reg) = curr_new_nt;
    end
    gene{1,mut_regions{r}} = {region_seq};
end

%Modifying phastCons and phyloP

if ~strcmp(mut_type,'snp')
    phastcons20 = gene.phastcons20{1};
    phastcons100 = gene.phastcons100{1};
    phylops20 = gene.phylops20{1};
    phylops100 = gene.phylops100{1};

    if strcmp(mut_type,'insert')
        phastcons20 = [phastcons20(1:min(mut_pos)-1); transpose(zeros(size(mut_pos))); phastcons20(min(mut_pos):end)];
        phastcons100 = [phastcons100(1:min(mut_pos)-1); transpose(zeros(size(mut_pos))); phastcons100(min(mut_pos):end)];
        phylops20 = [phylops20(1:min(mut_pos)-1); transpose(zeros(size(mut_pos))); phylops20(min(mut_pos):end)];
        phylops100 = [phylops100(1:min(mut_pos)-1); transpose(zeros(size(mut_pos))); phylops100(min(mut_pos):end)];
    elseif strcmp(mut_type,'del')
        phastcons20(mut_pos) = [];
        phastcons100(mut_pos) = [];
        phylops20(mut_pos) = [];
        phylops100(mut_pos) = [];
    end

    gene.phastcons20{1} = phastcons20;
    gene.phastcons100{1} = phastcons100;
    gene.phylops20{1} = phylops20;
    gene.phylops100{1} = phylops100;
end

gene.mut_ind(1) = mut_ind;
mut_gene = gene;

end