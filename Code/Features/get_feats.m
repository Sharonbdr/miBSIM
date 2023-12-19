function features = get_feats(UTR5,ORF,UTR3,region,seed_length,RNA_start,miRNA,seed_type,user,...
    with_cons_orf,phastcons100,phastcons20,phylops100,phylops20)

% UTR5/ORF/UTR3        =   Sequence of the region
% region               =   Of the binding site
% seed_length          =   6/7/8
% RNA_start            =   beginning of the target site in the complete mRNA
% miRNA                =   Sequence of the miRNA
% seed_type            =   '6mer'/'7mer-m8'/'7mer-A1'/'8mer'
% user                 =   user string, for the purpose of creating a scratch dir
% with_cons_orf        =   0/1 index for including conservation and ORF features, i.e. if the
%                          gene is human and endogenous or non-human and engineered.
% with_features        =   0/1 index for including calculated features in the output or only the
%                          predicted repression
% with_repression      =   0/1 index for including predicted repression in the output or only the
% features
% phastcons and phylop =   Needed if with_cons_orf == 1; phastCons and phyloP score vectors for the gene

regions = {'ORF','UTR3'};
seeds = {'6mer','7mer-A1','7mer-m8','8mer'};
seed_lens = [6,7,7,8];

r = find(strcmp(regions,region));
if isempty(r)
    error('region should be ORF/UTR3');
end

s = find(strcmp(seeds,seed_type));
if isempty(s)
    error('seed_type should be 6mer/7mer-A1/7mer-m8/8mer');
elseif seed_length ~= seed_lens(s)
    warning('seed_length doesn''t match seed_type');
end

%Creating local ViennaRNA folder
glob_path=[pwd,'/Thermo/ViennaRNA'];
loc_path = sprintf('/scratch/%s/',user);

if exist(loc_path,'dir') ~= 7
    mkdir(loc_path);
end

vienna_dir = [loc_path,'ViennaRNA'];
if exist(vienna_dir,'dir') == 7
    try 
        rmdir(vienna_dir,'s');
    catch
        warning(sprintf("could not delete vienna scratch folder at %s, please delete manualy", vienna_dir));
    end
end

copyfile(glob_path,vienna_dir);


%Calculating features
features = calc_gen_features(UTR5,ORF,UTR3,region,seed_length,RNA_start,miRNA,seed_type,glob_path,vienna_dir,with_cons_orf,...
    phastcons100,phastcons20,phylops100,phylops20);

% Removes ViennaRNA scratch directory
try 
    rmdir(vienna_dir,'s');
catch
    warning(sprintf("could not delete vienna scratch folder at %s, please delete manualy", vienna_dir));
end


end

