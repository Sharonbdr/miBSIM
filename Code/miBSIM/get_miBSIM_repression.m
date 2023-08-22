function []=get_miBSIM_repression(dataset) 
% Adds neighboring site information to feature set and corrects nominal
% repression predicted by the Intagrative model using miBSIM.
%
% dataset    =  A folder in Data that contains a genes table and miRs
%               Map Container.
%

% Variables
load(sprintf('../../Data/%s/genes.mat',dataset),'genes');
load(sprintf('../../Output/Integrative_output_%s.mat',dataset),'features');

% Adds neighboring site features
features = prepData(features);
save(sprintf('../../Output/check_%s.mat',dataset),'features');

% Predicts repression by miBSIM
% **NOTE** to get full feature set and per binding site repression return
% second output and save it.
[predictions,~] = correct_pred(features,genes);

save(sprintf('../../Output/%s_miBSIM_pred.mat',dataset),'predictions');
end