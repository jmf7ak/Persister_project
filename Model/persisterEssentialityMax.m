function [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMax(pa14, U5, U24, T5, T24)
% persisterEssentiality - integrates both RNA-sequencing and metabolomics data
% with the PA14 model to generate condition specific models. This script is
% based off of run_integratedEssentiality.
% 
% INPUTS:
%   pa14 - base pa14 model (biomass needs to be set to 'PA14_Biomass' and
%           media needs to be set to LB)
%   U5 - untreated at 5-hours condition-specific model
%   U24 - untreated at 24-hours condition-specific model
%   T5 - treated at 5-hours condition-specific model
%   T24 - treated at 24-hours condition-specific model
% 
% OUTPUTS:
%   model_ess_pa14 - essential genes for the pa14 model
%   model_ess_U5 - essential genes for the untreated at 5-hours model
%   model_ess_U24 - essential genes for the untreated at 24-hours model
%   model_ess_T5 - essential genes for the treated at 5-hours model
%   model_ess_T24 - essential genes for the treated at 24-hours model
%   common_ess - essential genes common to all conditions
%   unique_treated - essential genes unique to the treated conditions
%   unique_untreated- essential genes unique to the untreated conditions
%   unique_U5 - essential genes unique to the untreated at 5-hours model
%   unique_U24 - essential genes unique to the untreated at 24-hours model
%   unique_T5 - essential genes unique to the treated at 5-hours model
%   unique_T24 - essential genes unique to the treated at 24-hours model

%
% Anna Blazier, 2018-05-15

% NOTE: before running the script, run this:
% clear all
% close all
% clc
% initCobraToolbox()
% changeCobraSolver('gurobi5');
% model = d_xls2model_JAB('/Users/Anna/Dropbox/pseudomonas_model/PA14recon1_v24_published.xlsx');
% model = changeMedia_SEED(model, 1,'');
% model = changeObjective(model,'PA14_Biomass',1);
% sol = optimizeCbModel(model)

%% identify genes present in the model

deleted_T24 = [];
deleted_T5 = [];
deleted_U24 = [];
deleted_U5 = [];

present_T24 = [];
present_T5 = [];
present_U24 = [];
present_U5 = [];

deleted_ind_T24 = [];
deleted_ind_T5 = [];
deleted_ind_U24 = [];
deleted_ind_U5 = [];

for i = 1:length(pa14.genes)
    if isempty(findstr(char(T24.genes(i)),'_deleted')) == 1
        present_T24 = [present_T24,T24.genes(i)];
    else
        deleted_T24 = [deleted_T24,T24.genes(i)];
        deleted_ind_T24 = [deleted_ind_T24,i];
    end
    
    if isempty(findstr(char(T5.genes(i)),'_deleted')) == 1
        present_T5 = [present_T5,T5.genes(i)];
    else
        deleted_T5 = [deleted_T5,T5.genes(i)];
        deleted_ind_T5 = [deleted_ind_T5,i];
    end
    
    if isempty(findstr(char(U24.genes(i)),'_deleted')) == 1
        present_U24 = [present_U24,U24.genes(i)];
    else
        deleted_U24 = [deleted_U24,U24.genes(i)];
        deleted_ind_U24 = [deleted_ind_U24,i];
    end
    
    if isempty(findstr(char(U5.genes(i)),'_deleted')) == 1
        present_U5 = [present_U5,U5.genes(i)];
    else
        deleted_U5 = [deleted_U5,U5.genes(i)];
        deleted_ind_U5 = [deleted_ind_U5,i];
    end
    
end

%% perform single gene deletions

grRateKO_T24 = zeros(1,length(present_T24));
grRateKO_T5 = zeros(1,length(present_T5));
grRateKO_U24 = zeros(1,length(present_U24));
grRateKO_U5 = zeros(1,length(present_U5));
grRateKO_pa14 = zeros(1,length(pa14.genes));

for i = 1:length(present_T24)
    T24_deleted = deleteModelGenes(T24,present_T24(i));
    T24_deleted_sol = optimizeCbModel(T24_deleted,'max');
    T24_deleted_objFlux = T24_deleted_sol.f;
    grRateKO_T24(i) = T24_deleted_objFlux;
end
    
for i = 1:length(present_T5)
    T5_deleted = deleteModelGenes(T5,present_T5(i));
    T5_deleted_sol = optimizeCbModel(T5_deleted,'max');
    T5_deleted_objFlux = T5_deleted_sol.f;
    grRateKO_T5(i) = T5_deleted_objFlux;
end

for i = 1:length(present_U24)
    U24_deleted = deleteModelGenes(U24,present_U24(i));
    U24_deleted_sol = optimizeCbModel(U24_deleted,'max');
    U24_deleted_objFlux = U24_deleted_sol.f;
    grRateKO_U24(i) = U24_deleted_objFlux;
end

for i = 1:length(present_U5)
    U5_deleted = deleteModelGenes(U5,present_U5(i));
    U5_deleted_sol = optimizeCbModel(U5_deleted,'max');
    U5_deleted_objFlux = U5_deleted_sol.f;
    grRateKO_U5(i) = U5_deleted_objFlux;
end

for i = 1:length(pa14.genes)
    pa14_deleted = deleteModelGenes(pa14,pa14.genes(i));
    pa14_deleted_sol = optimizeCbModel(pa14_deleted,'max');
    pa14_deleted_objFlux = pa14_deleted_sol.f;
    grRateKO_pa14(i) = pa14_deleted_objFlux;
end

% [grRatio_T24, grRateKO_T24, grRateWT_T24, delRxns_T24, hasEffect_T24] = singleGeneDeletion(T24, 'FBA', present_T24);
% [grRatio_T5, grRateKO_T5, grRateWT_T5, delRxns_T5, hasEffect_T5] = singleGeneDeletion(T5, 'FBA', present_T5);
% [grRatio_U24, grRateKO_U24, grRateWT_U24, delRxns_U24, hasEffect_U24] = singleGeneDeletion(U24, 'FBA', present_U24);
% [grRatio_U5, grRateKO_U5, grRateWT_U5, delRxns_U5, hasEffect_U5] = singleGeneDeletion(U5, 'FBA', present_U5);
% [grRatio_pa14, grRateKO_pa14, grRateWT_pa14, delRxns_pa14, hasEffect_pa14] = singleGeneDeletion(pa14, 'FBA');

essentialInd_T24 = find(grRateKO_T24 < 0.0001);
essentialInd_T5 = find(grRateKO_T5 < 0.0001);
essentialInd_U24 = find(grRateKO_U24 < 0.0001);
essentialInd_U5 = find(grRateKO_U5 < 0.0001);
essentialInd_pa14 = find(grRateKO_pa14 < 0.0001);

NaN_T24 = find(isnan(grRateKO_T24));
NaN_T5 = find(isnan(grRateKO_T5));
NaN_U24 = find(isnan(grRateKO_U24));
NaN_U5 = find(isnan(grRateKO_U5));
NaN_pa14 = find(isnan(grRateKO_pa14));

essIdx_T24 = [essentialInd_T24 NaN_T24];
essIdx_T5 = [essentialInd_T5 NaN_T5];
essIdx_U24 = [essentialInd_U24 NaN_U24];
essIdx_U5 = [essentialInd_U5 NaN_U5];
essIdx_pa14 = [essentialInd_pa14 NaN_pa14];

model_ess_T24 = present_T24(essIdx_T24)';
model_ess_T5 = present_T5(essIdx_T5)';
model_ess_U24 = present_U24(essIdx_U24)';
model_ess_U5 = present_U5(essIdx_U5)';
model_ess_pa14 = pa14.genes(essIdx_pa14);

%% essential gene sets

common_ess = intersect(intersect(intersect(model_ess_T24,model_ess_T5),model_ess_U24), model_ess_U5);
treated_ess = intersect(model_ess_T24, model_ess_T5);
untreated_ess = intersect(model_ess_U24, model_ess_U5);

unique_treated = setdiff(setdiff(treated_ess, model_ess_U24), model_ess_U5);
unique_untreated = setdiff(setdiff(untreated_ess, model_ess_T24), model_ess_T5);

unique_T24 = setdiff(setdiff(setdiff(model_ess_T24, model_ess_T5), model_ess_U24), model_ess_U5);
unique_T5 = setdiff(setdiff(setdiff(model_ess_T5, model_ess_T24), model_ess_U24), model_ess_U5);
unique_U24 = setdiff(setdiff(setdiff(model_ess_U24, model_ess_T5), model_ess_T24), model_ess_U5);
unique_U5 = setdiff(setdiff(setdiff(model_ess_U5, model_ess_T5), model_ess_U24), model_ess_T24);
