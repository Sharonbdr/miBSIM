function [features] = predict_Intmodel(features)
% Predicts site repression based on feature set using the Integrative
% model.
%
% features  =  A table containing potential binding sites and their features
%

% Variables
load('../../Data/Integrative_model.mat','best_model');
feats_num=7:158;

regions = {'ORF','UTR3'};
seeds = {'6mer','7mer-A1','7mer-m8','8mer'};
with_riboseq = 1;

% Eliminates NaN and inf values
if with_riboseq
    tab_array = table2array(features(:,feats_num));
    [r,~] = find(isnan(tab_array) | isinf(tab_array) | imag(tab_array));
    nan_rows = unique(r);
    features(nan_rows,:) = [];
    clear r nan_rows
end


% Predicts per site repression
for r = 1:2 %region
    for s = 1:4 %seed
        region = regions{r};
        seed = seeds{s};
        a=best_model{r,s}.coeffs;
        b=best_model{r,s}.intercept;
        %specific site+region  features
        spec_f = features(strcmp(features.region,region) & strcmp(features.seed_type,seed),feats_num);

        %predicition based on features and current model
        new_pred = table2array(spec_f) * a + b;
        if sum(isnan(new_pred)) ~= 0
            error('NaN pred')
        end
        new_pred = min(new_pred,0);

        features.S(strcmp(features.region,region) & strcmp(features.seed_type,seed)) = new_pred;
    end
end


end