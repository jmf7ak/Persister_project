function [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(model, diffEX_file, RNA_pval, RNA_fc)
% persisterIntegration - integrates both RNA-sequencing and metabolomics data
% with the PA14 model to generate condition specific models. This script is
% based off of runMADE_180222.m and run_MetabolomicsIntegration.m
% 
% INPUTS:
%   model - base pa14 model (biomass needs to be set to 'PA14_Biomass' and
%           media needs to be set to LB)
%   diffEx_file - location of the differential expression information in 
%           .xls format
%   RNA_pval - p-value threshold for MADE
%   RNA_fc - log2FC threshold for MADE
% 
% OUTPUTS:
%   U5 - untreated at 5-hours condition-specific model
%   U24 - untreated at 24-hours condition-specific model
%   T5 - treated at 5-hours condition-specific model
%   T24 - treated at 24-hours condition-specific model
%   untreated_metabs_5_int - the integrated metabolites for the untreated
%           5-hours condition
%   persister_metabs_5_int - the integrated metabolites for the persister
%           5-hours condition
%   untreated_metabs_24_int - the integrated metabolites for the untreated
%           24-hours condition
%   persister_metabs_24_int - the integrated metabolites for the persister
%           24-hours condition
%
% Anna Blazier, 2018-05-15

% NOTE: before running the script, run this:
% clear all
% close all
% clc
% initCobraToolbox()
% changeCobraSolver('gurobi5');
% model = d_xls2model_JAB('/Users/Anna/Dropbox/pseudomonas_model/PA14recon1_v24_published.xlsx');
% model = changeMedia_SEED_new(model, 1,'');
% model = changeObjective(model,'PA14_Biomass',1);
% sol = optimizeCbModel(model)
% diffEX_file = 'C:\Users\Anna\Dropbox\Lab\Projects\Persisters\RNA_Sequencing\BIT_RNA_Sequencing\Second-run\Analysis\DESeq_dds_groupbatch_all2metab_integration_180221.xls';


%% Expression data load/format

% load the dataset
[ndata,tdata,rdata] = xlsread(diffEX_file);

% pull out the locus tags / gene ids
locus_tag = rdata(2:end,2);

% pull out the fold changes for my different conditions
fold_change_treated = cell2mat([rdata(2:end,3),rdata(2:end,5)]);
fold_change_untreated = cell2mat([rdata(2:end,7),rdata(2:end,9)]);

% pull out the p-values for my different conditions
p_values_treated = cell2mat([rdata(2:end,4),rdata(2:end,6)]);
p_values_untreated = cell2mat([rdata(2:end,8),rdata(2:end,10)]);

% set up transition matrices defining what comparisons are being made and
% how they relate to the order of the fold changes and p-values
T_treated = [ 0 1 2;
              0 0 0;
              0 0 0];
          
T_untreated = [ 0 1 2;
                0 0 0;
                0 0 0];
            
%% log2FC threshold implementation

p_values_treated(abs(fold_change_treated) < RNA_fc) = 10;
p_values_untreated(abs(fold_change_untreated) < RNA_fc) = 10;

%% Initialize tiger toolbox

% initiate the TIGER toolbox
start_tiger()

% convert the COBRA PA14 model to the TIGER PA14 model
modelT=cobra_to_tiger(model);

% set the solver to cplex
set_solver('cplex')

% confirm that flux through biomass is the same as COBRA using the TIGER
% toolbox
test=fba(modelT)

%% Setting up MADE

% set MADE solver options to display the results in real-time and to set a
% maximum run time
set_solver_option('Display','on'); 
set_solver_option('MaxTime',1000); 

% run MADE using the TIGER PA14 models, the fold change and p-value
% information, and the specified transition matrix
made_sol_treated = made(modelT, fold_change_treated, p_values_treated, 'gene_names', locus_tag, 'transition_matrix',T_treated, 'log_fold_change', 'true', 'verify', 'true', 'p_thresh', RNA_pval, 'obj_frac', 0.1);
made_sol_untreated = made(modelT, fold_change_untreated, p_values_untreated, 'gene_names', locus_tag, 'transition_matrix',T_untreated, 'log_fold_change', 'true', 'verify', 'true', 'p_thresh', RNA_pval, 'obj_frac', 0.1);

% display the MADE results
show_made_results(made_sol_treated,'all');
show_made_results(made_sol_untreated,'all');

%% Remove transition genes - treated
% the final step is to take the gene states suggested by MADE and implement
% them in COBRA models by performing gene knockouts. For some reason, MADE
% will randomly turn genes off completely across the different conditions.
% We don't want to include these in our integration since they aren't based
% on the differential expression analysis. Instead, we only want to include
% the gene states for which there is a "transition" (i.e., off in one
% condition and on in another). 

% initiate the COBRA toolbox
initCobraToolbox()
changeCobraSolver('gurobi');

% for each gene, sum the gene state activity across the conditions, this
% will be used to determine if a transition is happening or not
rowSum = [];
rowSum = sum(made_sol_treated.gene_states');

% initiate empty variables
transitionIndex_treated = [];
transitionGene_treated = [];
transitionState_treated = [];

% determine whether a gene undergoes a transition. this will need to be
% modified based on how many conditions are being considered.
for i = 1:length(rowSum)
    if rowSum(i) ~= 0 && rowSum(i) ~= 3
        transitionIndex_treated = [transitionIndex_treated,i];
        transitionGene_treated = [transitionGene_treated,made_sol_treated.genes(i)];
        transitionState_treated = [transitionState_treated; made_sol_treated.gene_states(i,1:3)];
    end
end

size_treated = size(transitionState_treated);

% initiate empty variables
treated_0_transition_Gene = [];
treated_5_transition_Gene = [];
treated_24_transition_Gene = [];

% get a list of the genes that need to be turned off for each condition
% based on whether or not a transition is happening
for i = 1:size_treated(1)
    if transitionState_treated(i,1) == 0
        treated_0_transition_Gene = [treated_0_transition_Gene,transitionGene_treated(i)];
    end
    if transitionState_treated(i,2) == 0
        treated_5_transition_Gene = [treated_5_transition_Gene,transitionGene_treated(i)];
    end
    if transitionState_treated(i,3) == 0
        treated_24_transition_Gene = [treated_24_transition_Gene,transitionGene_treated(i)];
    end
end

% create condition-specific COBRA models by deleting genes from the PA14
% model based on the transition gene state information. 

if isempty(treated_0_transition_Gene) ~= 1
    T0 = deleteModelGenes(model,treated_0_transition_Gene);
else
    T0 = model;
end
T0_transition_sol=optimizeCbModel(T0);

if isempty(treated_5_transition_Gene) ~= 1
    T5 = deleteModelGenes(model,treated_5_transition_Gene);
else
    T5 = model;
end
T5_transition_sol=optimizeCbModel(T5);

if isempty(treated_24_transition_Gene) ~= 1
    T24 = deleteModelGenes(model,treated_24_transition_Gene);
else
    T24 = model;
end
T24_transition_sol=optimizeCbModel(T24);

%% Remove transition genes - untreated

rowSum = [];
rowSum = sum(made_sol_untreated.gene_states');

transitionIndex_untreated = [];
transitionGene_untreated = [];
transitionState_untreated = [];

for i = 1:length(rowSum)
    if rowSum(i) ~= 0 && rowSum(i) ~= 3
        transitionIndex_untreated = [transitionIndex_untreated,i];
        transitionGene_untreated = [transitionGene_untreated,made_sol_untreated.genes(i)];
        transitionState_untreated = [transitionState_untreated; made_sol_untreated.gene_states(i,1:3)];
    end
end

size_untreated = size(transitionState_untreated);

untreated_0_transition_Gene = [];
untreated_5_transition_Gene = [];
untreated_24_transition_Gene = [];

for i = 1:size_untreated(1)
    if transitionState_untreated(i,1) == 0
        untreated_0_transition_Gene = [untreated_0_transition_Gene,transitionGene_untreated(i)];
    end
    if transitionState_untreated(i,2) == 0
        untreated_5_transition_Gene = [untreated_5_transition_Gene,transitionGene_untreated(i)];
    end
    if transitionState_untreated(i,3) == 0
        untreated_24_transition_Gene = [untreated_24_transition_Gene,transitionGene_untreated(i)];
    end
end

if isempty(untreated_0_transition_Gene) ~= 1
    U0 = deleteModelGenes(model,untreated_0_transition_Gene);
else
    U0 = model;
end
U0_transition_sol=optimizeCbModel(U0);

if isempty(untreated_5_transition_Gene) ~= 1
    U5 = deleteModelGenes(model,untreated_5_transition_Gene);
else
    U5 = model;
end
U5_transition_sol=optimizeCbModel(U5);

if isempty(untreated_24_transition_Gene) ~= 1
    U24 = deleteModelGenes(model,untreated_24_transition_Gene);
else
    U24 = model;
end
U24_transition_sol=optimizeCbModel(U24);

%% Metabolites for integration
% untreated metabolites
untreated_metabs_5 = {
    'cpd02095[cytosol]'; % 2,4-diaminobutyrate %%% flux in both U5 and pa14
    'cpd00117[extracellular]'; % D-alanine %%% NO FLUX in U5 or pa14
    'cpd00035[extracellular]'; % L-alanine %%% flux in both U5 and pa14
    'cpd11583[extracellular]'; % Ala-Leu %%% NO FLUX in U5 or pa14
    'cpd00024[extracellular]'; % 2-Oxoglutarate %%% flux in both U5 and pa14
    'cpd00132[extracellular]'; % asparagine %%% flux in both U5 and pa14
    'cpd00041[extracellular]'; % aspartate %%% flux in both U5 and pa14
    'cpd00098[extracellular]'; % choline %%% NO FLUX in U5 or pa14
    'cpd00137[extracellular]'; % citrate %%% flux in both U5 and pa14
    'cpd00053[extracellular]'; % glutamine %%% NO FLUX in U5 or pa14
    'cpd15604[extracellular]'; % Gly-Leu %%% NO FLUX in U5 or pa14
    'cpd11584[extracellular]'; % Ala-His %%% NO FLUX in U5 or pa14
    'cpd00790[cytosol]'; % O-Acetyl-L-homoserine %%% flux in both U5 and pa14
    'cpd15605[extracellular]'; % Gly-Phe %%% NO FLUX in U5 or pa14
    'cpd00129[extracellular]'; % proline %%% flux in both U5 and pa14
    'cpd00599[cytosol]'; % salicylate %%% flux in both U5 and pa14
    'cpd00523[cytosol]'; % trehalose-6-phosphate %%% NO FLUX in U5 or pa14
    'cpd15606[extracellular]'; % Gly-Tyr %%% NO FLUX in U5 or pa14
    };

% persister metabolites
persister_metabs_5 = {
    'cpd00169[cytosol]'; % 3-phosphoglycerate %%% flux in both T5 and pa14
    'cpd00868[cytosol]'; % p-hydroxyphenylpyruvate %%% flux in both T5 and pa14
    'cpd00331[extracellular]'; % aconitate %%% flux in both T5 and pa14
    'cpd00018[cytosol]'; % AMP %%% flux in both T5 and pa14
    'cpd00051[extracellular]'; % arginine %%% NO FLUX in T5 or pa14
    'cpd01155[cytosol]'; % cadaverine %%% flux in pa14 BUT NOT T5
    'cpd00581[extracellular]'; % Urocanate %%% flux in both T5 and pa14
    'cpd00106[extracellular]'; % fumarate %%% flux in both T5 and pa14
    'cpd00379[extracellular]'; % glutarate %%% flux in both T5 and pa14
    'cpd00080[extracellular]'; % glycerol-3-phosphate %%% flux in pa14 BUT NOT T5
    'cpd00873[cytosol]'; % nicotinamide ribonucleotide %%% flux in pa14 BUT NOT T5
    'cpd00247[extracellular]'; % orotate %%% flux in pa14 BUT NOT T5
    'cpd00061[cytosol]'; % phosphoenolpyruvate %%% flux in both T5 and pa14
    'cpd00036[extracellular]'; % succinate %%% flux in both T5 and pa14
    'cpd00092[extracellular]'; % uracil %%% flux in both T5 and pa14
    'cpd00309[extracellular]'; % XAN %%% flux in pa14 but NOT T5
    };

untreated_metabs_24 = {
    'cpd02095[cytosol]'; % 2,4-diaminobutyrate %%% flux in both U24 and pa14
    'cpd00876[cytosol]'; % 3-hydroxyisobutyrate %%% flux in pa14 BUT NOT U24
    'cpd00169[cytosol]'; % 3-phosphoglycerate %%% flux in both U24 and pa14
    'cpd00136[extracellular]'; % 4-hydroxybenzoate %%% flux in both U24 and pa14
    'cpd00489[extracellular]'; % 4-hydroxyphenylacetate %%% NO FLUX in U24 or pa14
    'cpd01293[extracellular]'; % 5-oxoproline %%% NO FLUX in U24 or pa14
    'cpd01155[cytosol]'; % cadaverine %%% flux in both U24 and pa14
    'cpd00015[cytosol]'; % FAD %%% flux in pa14 BUT NOT U24
    'cpd00027[extracellular]'; % glucose %%% NO FLUX in U24 or pa14
    'cpd00053[extracellular]'; % glutamine %%% NO FLUX in U24 or pa14
    'cpd00033[extracellular]'; % glycine %%% flux in both U24 and pa14
    'cpd00126[cytosol]'; % GMP %%% flux in both U24 and pa14
    'cpd00119[extracellular]'; % histidine %%% flux in pa14 but not U24
    'cpd00227[cytosol]'; % homoserine %%% flux in pa14 BUT NOT U24
    'cpd00322[extracellular]'; % isoleucine %%% flux in pa14 but not U24
    'cpd00107[extracellular]'; % leucine %%% flux in pa14 but not U24
    'cpd00039[extracellular]'; % lysine %%% flux in both U24 and pa14
    'cpd00342[cytosol]'; % N-Acetylornithine %%% flux in both U24 and pa14
    'cpd00790[cytosol]'; % O-Acetyl-L-homoserine %%% flux in both U24 and pa14
    'cpd00064[extracellular]'; % ornithine %%% flux in both U24 and pa14
    'cpd00247[extracellular]'; % orotate %%% flux in pa14 BUT NOT U24
    'cpd02333[cytosol]'; % quinolate %%% flux in both U24 and pa14
    'cpd00105[extracellular]'; % ribose %%% flux in pa14 but not U24
    'cpd00550[extracellular]'; % D-serine %%% NO FLUX in U24 or pa14
    'cpd00054[extracellular]'; % L-serine %%% flux in both U24 and pa14
    'cpd00048[extracellular]'; % sulfate %%% NO FLUX in U24 or pa14
    'cpd00161[extracellular]'; % threonine %%% flux in pa14 but not U24
    'cpd00184[extracellular]'; % thymidine %%% flux in both U24 and pa14
    'cpd00073[extracellular]'; % urea %%% flux in both U24 and pa14
    };

persister_metabs_24 = {
    'cpd00738[cytosol]'; % phosphoserine %%% flux in both T24 and pa14
    'cpd00147[extracellular]'; % 5-methylthioadenosine %%% NO FLUX in T24 or pa14
    'cpd11583[extracellular]'; % Ala-Leu %%% NO FLUX in T24 or pa14
    'cpd00024[extracellular]'; % 2-Oxoglutarate %%% flux in both T24 and pa14
    'cpd00051[extracellular]'; % arginine %%% NO FLUX in T24 or pa14
    'cpd00041[extracellular]'; % aspartate %%% flux in pa14 but not T24
    'cpd00581[extracellular]'; % urocanate %%% flux in both T24 and pa14
    'cpd00046[cytosol]'; % CMP %%% flux in both T24 and pa14
    'cpd00080[extracellular]'; % glycerol-3-phosphate %%% flux in both T24 and pa14
    'cpd00061[cytosol]'; % phosphoenolpyruvate %%% flux in both T24 and pa14
    'cpd00129[extracellular]'; % proline %%% flux in pa14 but not T24
    'cpd00020[extracellular]'; % pyruvate %%% flux in both T24 and pa14
    'cpd00794[extracellular]'; % TRHL %%% NO FLUX in T24 or pa14
    'cpd00309[extracellular]'; % XAN %%% flux in both T24 and pa14
    };

%% Force production through single metabolite to test model

rxnIdx = zeros(length(untreated_metabs_5),1);
untreated_metabs_5_int = {};

for i = 1:length(untreated_metabs_5)
    
    % create test models to see if it's even possible to integrate the
    % metab
    U5_test = U5;
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(untreated_metabs_5{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(U5.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        U5_test = addDemandReaction(U5_test,untreated_metabs_5{i});
        dmRxnName = strcat('DM_',untreated_metabs_5{i});
        dmRxnIdx = find(ismember(U5_test.rxns,dmRxnName));
        U5_test.lb(dmRxnIdx) = 0.0001;
        U5_test.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        U5_test.lb(exRxnIdx) = 0.0001;
        U5_test.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
    
    % measure flux through biomass when forcing production of the single
    % metabolite
    sol = optimizeCbModel(U5_test);
    if sol.f ~= 0 && ~isnan(sol.f)
        untreated_metabs_5_int{end+1} = untreated_metabs_5{i};
    end
    
end

rxnIdx = zeros(length(persister_metabs_5),1);
persister_metabs_5_int = {};

for i = 1:length(persister_metabs_5)
    
    % create test models to see if it's even possible to integrate the
    % metab
    T5_test = T5;
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(persister_metabs_5{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(T5.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        T5_test = addDemandReaction(T5_test,persister_metabs_5{i});
        dmRxnName = strcat('DM_',persister_metabs_5{i});
        dmRxnIdx = find(ismember(T5_test.rxns,dmRxnName));
        T5_test.lb(dmRxnIdx) = 0.0001;
        T5_test.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        T5_test.lb(exRxnIdx) = 0.0001;
        T5_test.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
    
    % measure flux through biomass when forcing production of the single
    % metabolite
    sol = optimizeCbModel(T5_test);
    if sol.f ~= 0 && ~isnan(sol.f)
        persister_metabs_5_int{end+1} = persister_metabs_5{i};
    end
    
end    
    
rxnIdx = zeros(length(untreated_metabs_24),1);
untreated_metabs_24_int = {};

for i = 1:length(untreated_metabs_24)
    
    % create test models to see if it's even possible to integrate the
    % metab
    U24_test = U24;
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(untreated_metabs_24{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(U24.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        U24_test = addDemandReaction(U24_test,untreated_metabs_24{i});
        dmRxnName = strcat('DM_',untreated_metabs_24{i});
        dmRxnIdx = find(ismember(U24_test.rxns,dmRxnName));
        U24_test.lb(dmRxnIdx) = 0.0001;
        U24_test.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        U24_test.lb(exRxnIdx) = 0.0001;
        U24_test.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
    
    % measure flux through biomass when forcing production of the single
    % metabolite
    sol = optimizeCbModel(U24_test);
    if sol.f ~= 0 && ~isnan(sol.f)
        untreated_metabs_24_int{end+1} = untreated_metabs_24{i};
    end
    
end

rxnIdx = zeros(length(persister_metabs_24),1);
persister_metabs_24_int = {};

for i = 1:length(persister_metabs_24)
    
    % create test models to see if it's even possible to integrate the
    % metab
    T24_test = T24;
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(persister_metabs_24{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(T24.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        T24_test = addDemandReaction(T24_test,persister_metabs_24{i});
        dmRxnName = strcat('DM_',persister_metabs_24{i});
        dmRxnIdx = find(ismember(T24_test.rxns,dmRxnName));
        T24_test.lb(dmRxnIdx) = 0.0001;
        T24_test.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        T24_test.lb(exRxnIdx) = 0.0001;
        T24_test.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
    
    % measure flux through biomass when forcing production of the single
    % metabolite
    sol = optimizeCbModel(T24_test);
    if sol.f ~= 0 && ~isnan(sol.f)
        persister_metabs_24_int{end+1} = persister_metabs_24{i};
    end
    
end     

%% force production of metabolites

% untreated at 5 hours
rxnIdx = zeros(length(untreated_metabs_5_int),1);

for i = 1:length(untreated_metabs_5_int)
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(untreated_metabs_5_int{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(U5.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        U5 = addDemandReaction(U5,untreated_metabs_5_int{i});
        dmRxnName = strcat('DM_',untreated_metabs_5_int{i});
        dmRxnIdx = find(ismember(U5.rxns,dmRxnName));
        U5.lb(dmRxnIdx) = 0.0001;
        U5.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        U5.lb(exRxnIdx) = 0.0001;
        U5.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
end

% persister at 5 hours
rxnIdx = zeros(length(persister_metabs_5_int),1);

for i = 1:length(persister_metabs_5_int)
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(persister_metabs_5_int{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(T5.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        T5 = addDemandReaction(T5,persister_metabs_5_int{i});
        dmRxnName = strcat('DM_',persister_metabs_5_int{i});
        dmRxnIdx = find(ismember(T5.rxns,dmRxnName));
        T5.lb(dmRxnIdx) = 0.0001;
        T5.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        T5.lb(exRxnIdx) = 0.0001;
        T5.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
end

% untreated at 24 hours
rxnIdx = zeros(length(untreated_metabs_24_int),1);

for i = 1:length(untreated_metabs_24_int)
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(untreated_metabs_24_int{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(U24.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        U24 = addDemandReaction(U24,untreated_metabs_24_int{i});
        dmRxnName = strcat('DM_',untreated_metabs_24_int{i});
        dmRxnIdx = find(ismember(U24.rxns,dmRxnName));
        U24.lb(dmRxnIdx) = 0.0001;
        U24.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        U24.lb(exRxnIdx) = 0.0001;
        U24.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
end

% persister at 24 hours
rxnIdx = zeros(length(persister_metabs_24_int),1);

for i = 1:length(persister_metabs_24_int)
    
    % check to see if an exchange reaction for the metab exists
    split = strsplit(persister_metabs_24_int{i},'[');
    exRxnName = strcat('EX_',split(1),'_e');
    exRxnIdx = find(ismember(T24.rxns,exRxnName));
    
    % if there is no exchange reaction, create a demand reaction
    if isempty(exRxnIdx) == 1
        T24 = addDemandReaction(T24,persister_metabs_24_int{i});
        dmRxnName = strcat('DM_',persister_metabs_24_int{i});
        dmRxnIdx = find(ismember(T24.rxns,dmRxnName));
        T24.lb(dmRxnIdx) = 0.0001;
        T24.ub(dmRxnIdx) = 1000;
        rxnIdx(i) = dmRxnIdx;
        
    % if there is an exchange reaction, modify its bounds to force production    
    else
        T24.lb(exRxnIdx) = 0.0001;
        T24.ub(exRxnIdx) = 1000;
        rxnIdx(i) = exRxnIdx;
    end
end
