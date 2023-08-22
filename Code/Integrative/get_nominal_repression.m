function [] = get_nominal_repression(dataset)
% Gathers base feature set and biochemical features to one feature table,
% and predicts per site repression using the Integrative model.
%
% dataset    =  A folder in Data that contains a genes table and miRs
%               Map Container.
%

% Variables
load(sprintf('../../Data/%s/genes.mat',dataset),'genes');
load(sprintf('../../Output/Output_%s.mat',dataset),'Output');
cnn_results_path=sprintf('../Features/Biochemical/reports/outputs/predictions/%s/',dataset);

% Combine all features
features = add_bio_features(genes,Output,cnn_results_path);

% predict Integrative model
features = predict_Intmodel(features);

save(sprintf('../../Output/Integrative_output_%s.mat',dataset),'features');
end