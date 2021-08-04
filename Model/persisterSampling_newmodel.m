%% initialize environment

clear all
close all
clc
initCobraToolbox()
changeCobraSolver('gurobi');
%%
% Windows:
pa14 = readCbModel('iPau21');
pa14 = creategrRulesField(pa14);
pa14 = changeMedia_SEED_new(pa14, 1,'');
pa14 = changeObjective(pa14,'PA14_Biomass',1);
sol = optimizeCbModel(pa14);

%%
% Windows:
diffEX_file = 'C:\Users\jmfic\OneDrive\Documents\21_SUMMA\Persister\model\data\DESeq_dds_groupbatch_all2metab_integration_joe.csv';


%% Create basal biomass models

pa14Bio = changeObjective(pa14,'PA14_Biomass',1);
pa14Bio_sol = optimizeCbModel(pa14Bio);
pa14Bio = creategrRulesField(pa14Bio);

RNA_pval = [0.5];
RNA_fc = [0];

[U5Bio, U24Bio, T5Bio, T24Bio, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration_newmodel(pa14Bio, diffEX_file, RNA_pval, RNA_fc);

U5Bio_sol = optimizeCbModel(U5Bio);
T5Bio_sol = optimizeCbModel(T5Bio);
U24Bio_sol = optimizeCbModel(U24Bio);
T24Bio_sol = optimizeCbModel(T24Bio);
        
U5Bio.ub(1487) = 0.1*U5Bio_sol.f;
T5Bio.ub(1487) = 0.1*T5Bio_sol.f;
U24Bio.ub(1487) = 0.1*U24Bio_sol.f;
T24Bio.ub(1487) = 0.1*T24Bio_sol.f;
pa14Bio.ub(1487) = 0.1*pa14Bio_sol.f;

%% Flux sampling on basal biomass models
% Note, in order to run this next chunk, I needed to change the directory 
% to optGpSampler_1.1_Matlab

flux_U5Bio = optGpSampler(U5Bio,[],3000,1,1,'gurobi',0);
flux_U24Bio = optGpSampler(U24Bio,[],3000,1,1,'gurobi',0);
flux_T5Bio = optGpSampler(T5Bio,[],3000,1,1,'gurobi',0);
flux_T24Bio = optGpSampler(T24Bio,[],3000,1,1,'gurobi',0);

fluxDist_U5Bio = zeros(length(U5Bio.rxns),3000);
fluxMean_U5Bio = zeros(length(U5Bio.rxns),1);
fluxMedian_U5Bio = zeros(length(U5Bio.rxns),1);
for i = 1:length(U5Bio.rxns)
    fluxDist_U5Bio(i,:) = flux_U5Bio.points(findRxnIDs(U5Bio,U5Bio.rxns(i)),:);
    fluxMean_U5Bio(i) = mean(fluxDist_U5Bio(i,:));
    fluxMedian_U5Bio(i) = median(fluxDist_U5Bio(i,:));
end

fluxDist_U24Bio = zeros(length(U24Bio.rxns),3000);
fluxMean_U24Bio = zeros(length(U24Bio.rxns),1);
fluxMedian_U24Bio = zeros(length(U24Bio.rxns),1);
for i = 1:length(U24Bio.rxns)
    fluxDist_U24Bio(i,:) = flux_U24Bio.points(findRxnIDs(U24Bio,U24Bio.rxns(i)),:);
    fluxMean_U24Bio(i) = mean(fluxDist_U24Bio(i,:));
    fluxMedian_U24Bio(i) = median(fluxDist_U24Bio(i,:));
end

fluxDist_T5Bio = zeros(length(T5Bio.rxns),3000);
fluxMean_T5Bio = zeros(length(T5Bio.rxns),1);
fluxMedian_T5Bio = zeros(length(T5Bio.rxns),1);
for i = 1:length(T5Bio.rxns)
    fluxDist_T5Bio(i,:) = flux_T5Bio.points(findRxnIDs(T5Bio,T5Bio.rxns(i)),:);
    fluxMean_T5Bio(i) = mean(fluxDist_T5Bio(i,:));
    fluxMedian_T5Bio(i) = median(fluxDist_T5Bio(i,:));
end

fluxDist_T24Bio = zeros(length(T24Bio.rxns),3000);
fluxMean_T24Bio = zeros(length(T24Bio.rxns),1);
fluxMedian_T24Bio = zeros(length(T24Bio.rxns),1);
for i = 1:length(T24Bio.rxns)
    fluxDist_T24Bio(i,:) = flux_T24Bio.points(findRxnIDs(T24Bio,T24Bio.rxns(i)),:);
    fluxMean_T24Bio(i) = mean(fluxDist_T24Bio(i,:));
    fluxMedian_T24Bio(i) = median(fluxDist_T24Bio(i,:));
end

% for each reaction, find the proportion of fluxes greater than 0.00001
proportion = zeros(length(pa14.rxns),4);

for i = 1:length(pa14.rxns)
    activeSamples_U5Bio = find(fluxDist_U5Bio(i,:) > 0.0001);
    activeSamples_U24Bio = find(fluxDist_U24Bio(i,:) > 0.0001);
    activeSamples_T5Bio = find(fluxDist_T5Bio(i,:) > 0.0001);
    activeSamples_T24Bio = find(fluxDist_T24Bio(i,:) > 0.0001);

    proportion(i,1) = length(activeSamples_U5Bio)/length(fluxDist_U5Bio);
    proportion(i,2) = length(activeSamples_U24Bio)/length(fluxDist_U24Bio);
    proportion(i,3) = length(activeSamples_T5Bio)/length(fluxDist_T5Bio);
    proportion(i,4) = length(activeSamples_T24Bio)/length(fluxDist_T24Bio);
end

% tidy data for R
% NOTE THAT FOR THIS ANALYSIS, WE ARE ONLY CONSIDERING REACTIONS THAT ARE
% IN ALL OF THE MODELS...SO THE FORCED METABOLITE PRODUCTION REACTIONS WILL
% NOT BE CONSIDERED...PERHAPS WANT TO CHANGE THIS FOR THE CORRELATION
% STUFF?

fluxDist_U5Bio_R = [];
fluxDist_U24Bio_R = [];
fluxDist_T5Bio_R = [];
fluxDist_T24Bio_R = [];
fluxDistBio = [];

for i = 1:length(pa14.rxns)
    fluxDist_U5Bio_R = [fluxDist_U5Bio(i,:)'];
    fluxDist_U24Bio_R = [fluxDist_U24Bio(i,:)'];
    fluxDist_T5Bio_R = [fluxDist_T5Bio(i,:)'];
    fluxDist_T24Bio_R = [fluxDist_T24Bio(i,:)'];
    fluxDistBio = [fluxDistBio;
                fluxDist_U5Bio_R;
                fluxDist_U24Bio_R;
                fluxDist_T5Bio_R;
                fluxDist_T24Bio_R];
end

fluxDist_rxns_R = [];
fluxDist_conditions_R = [];

conditions_R = [repmat({'U5'},[3000 1]);
                repmat({'U24'},[3000 1]);
                repmat({'T5'}, [3000 1]);
                repmat({'T24'}, [3000 1])];

for i = 1:length(pa14.rxns)
    rxns_mat = repmat(pa14.rxns(i),(3000*4),1);
    if i == 1
        fluxDist_rxns_R = rxns_mat;
    else
        fluxDist_rxns_R = [fluxDist_rxns_R;rxns_mat];
    end
    
    fluxDist_conditions_R = [fluxDist_conditions_R; conditions_R];
end

% I can't save these vectors in excel cause they are too large.
% The max length is 1048576

fluxDist_rxns_R_1 = fluxDist_rxns_R(1:1048576);
fluxDist_rxns_R_2 = fluxDist_rxns_R(1048577:2097152);
fluxDist_rxns_R_3 = fluxDist_rxns_R(2097153:3145728);
fluxDist_rxns_R_4 = fluxDist_rxns_R(3145729:4194300);
fluxDist_rxns_R_5 = fluxDist_rxns_R(4194301:5242871);
fluxDist_rxns_R_6 = fluxDist_rxns_R(5242872:6291447);
fluxDist_rxns_R_7 = fluxDist_rxns_R(6291448:7340023);
fluxDist_rxns_R_8 = fluxDist_rxns_R(7340024:8388509);
fluxDist_rxns_R_9 = fluxDist_rxns_R(8388510:9437085);
fluxDist_rxns_R_10 = fluxDist_rxns_R(9437086:10485661);
fluxDist_rxns_R_11 = fluxDist_rxns_R(10485662:11534237);
fluxDist_rxns_R_12 = fluxDist_rxns_R(11534238:12582807);
fluxDist_rxns_R_13 = fluxDist_rxns_R(12582808:13631383);
fluxDist_rxns_R_14 = fluxDist_rxns_R(13631384:14679959);
fluxDist_rxns_R_15 = fluxDist_rxns_R(14679960:15728535);
fluxDist_rxns_R_16 = fluxDist_rxns_R(15728536:16777111);
fluxDist_rxns_R_17 = fluxDist_rxns_R(16777112:17825687);
fluxDist_rxns_R_18 = fluxDist_rxns_R(17825688:end);

fluxDist_conditions_R_1 = fluxDist_conditions_R(1:1048576);
fluxDist_conditions_R_2 = fluxDist_conditions_R(1048577:2097152);
fluxDist_conditions_R_3 = fluxDist_conditions_R(2097153:3145728);
fluxDist_conditions_R_4 = fluxDist_conditions_R(3145729:4194300);
fluxDist_conditions_R_5 = fluxDist_conditions_R(4194301:5242871);
fluxDist_conditions_R_6 = fluxDist_conditions_R(5242872:6291447);
fluxDist_conditions_R_7 = fluxDist_conditions_R(6291448:7340023);
fluxDist_conditions_R_8 = fluxDist_conditions_R(7340024:8388509);
fluxDist_conditions_R_9 = fluxDist_conditions_R(8388510:9437085);
fluxDist_conditions_R_10 = fluxDist_conditions_R(9437086:10485661);
fluxDist_conditions_R_11 = fluxDist_conditions_R(10485662:11534237);
fluxDist_conditions_R_12 = fluxDist_conditions_R(11534238:12582807);
fluxDist_conditions_R_13 = fluxDist_conditions_R(12582808:13631383);
fluxDist_conditions_R_14 = fluxDist_conditions_R(13631384:14679959);
fluxDist_conditions_R_15 = fluxDist_conditions_R(14679960:15728535);
fluxDist_conditions_R_16 = fluxDist_conditions_R(15728536:16777111);
fluxDist_conditions_R_17 = fluxDist_conditions_R(16777112:17825687);
fluxDist_conditions_R_18 = fluxDist_conditions_R(17825688:end);

fluxDist_R_1 = fluxDistBio(1:1048576);
fluxDist_R_2 = fluxDistBio(1048577:2097152);
fluxDist_R_3 = fluxDistBio(2097153:3145728);
fluxDist_R_4 = fluxDistBio(3145729:4194300);
fluxDist_R_5 = fluxDistBio(4194301:5242871);
fluxDist_R_6 = fluxDistBio(5242872:6291447);
fluxDist_R_7 = fluxDistBio(6291448:7340023);
fluxDist_R_8 = fluxDistBio(7340024:8388509);
fluxDist_R_9 = fluxDistBio(8388510:9437085);
fluxDist_R_10 = fluxDistBio(9437086:10485661);
fluxDist_R_11 = fluxDistBio(10485662:11534237);
fluxDist_R_12 = fluxDistBio(11534238:12582807);
fluxDist_R_13 = fluxDistBio(12582808:13631383);
fluxDist_R_14 = fluxDistBio(13631384:14679959);
fluxDist_R_15 = fluxDistBio(14679960:15728535);
fluxDist_R_16 = fluxDistBio(15728536:16777111);
fluxDist_R_17 = fluxDistBio(16777112:17825687);
fluxDist_R_18 = fluxDistBio(17825688:end);

xlswrite('Bio_sampling1.xlsx',fluxDist_rxns_R_1,'Sheet1','A1')
xlswrite('Bio_sampling1.xlsx',fluxDist_conditions_R_1,'Sheet1','B1')
xlswrite('Bio_sampling1.xlsx',fluxDist_R_1,'Sheet1','C1')

xlswrite('Bio_sampling2.xlsx',fluxDist_rxns_R_2,'Sheet1','A1')
xlswrite('Bio_sampling2.xlsx',fluxDist_conditions_R_2,'Sheet1','B1')
xlswrite('Bio_sampling2.xlsx',fluxDist_R_2,'Sheet1','C1')

xlswrite('Bio_sampling3.xlsx',fluxDist_rxns_R_3,'Sheet1','A1')
xlswrite('Bio_sampling3.xlsx',fluxDist_conditions_R_3,'Sheet1','B1')
xlswrite('Bio_sampling3.xlsx',fluxDist_R_3,'Sheet1','C1')

xlswrite('Bio_sampling4.xlsx',fluxDist_rxns_R_4,'Sheet1','A1')
xlswrite('Bio_sampling4.xlsx',fluxDist_conditions_R_4,'Sheet1','B1')
xlswrite('Bio_sampling4.xlsx',fluxDist_R_4,'Sheet1','C1')

xlswrite('Bio_sampling5.xlsx',fluxDist_rxns_R_5,'Sheet1','A1')
xlswrite('Bio_sampling5.xlsx',fluxDist_conditions_R_5,'Sheet1','B1')
xlswrite('Bio_sampling5.xlsx',fluxDist_R_5,'Sheet1','C1')

xlswrite('Bio_sampling6.xlsx',fluxDist_rxns_R_6,'Sheet1','A1')
xlswrite('Bio_sampling6.xlsx',fluxDist_conditions_R_6,'Sheet1','B1')
xlswrite('Bio_sampling6.xlsx',fluxDist_R_6,'Sheet1','C1')

xlswrite('Bio_sampling7.xlsx',fluxDist_rxns_R_7,'Sheet1','A1')
xlswrite('Bio_sampling7.xlsx',fluxDist_conditions_R_7,'Sheet1','B1')
xlswrite('Bio_sampling7.xlsx',fluxDist_R_7,'Sheet1','C1')

xlswrite('Bio_sampling8.xlsx',fluxDist_rxns_R_8,'Sheet1','A1')
xlswrite('Bio_sampling8.xlsx',fluxDist_conditions_R_8,'Sheet1','B1')
xlswrite('Bio_sampling8.xlsx',fluxDist_R_8,'Sheet1','C1')

xlswrite('Bio_sampling9.xlsx',fluxDist_rxns_R_9,'Sheet1','A1')
xlswrite('Bio_sampling9.xlsx',fluxDist_conditions_R_9,'Sheet1','B1')
xlswrite('Bio_sampling9.xlsx',fluxDist_R_9,'Sheet1','C1')

xlswrite('Bio_sampling10.xlsx',fluxDist_rxns_R_10,'Sheet1','A1')
xlswrite('Bio_sampling10.xlsx',fluxDist_conditions_R_10,'Sheet1','B1')
xlswrite('Bio_sampling10.xlsx',fluxDist_R_10,'Sheet1','C1')

xlswrite('Bio_sampling11.xlsx',fluxDist_rxns_R_11,'Sheet1','A1')
xlswrite('Bio_sampling11.xlsx',fluxDist_conditions_R_11,'Sheet1','B1')
xlswrite('Bio_sampling11.xlsx',fluxDist_R_11,'Sheet1','C1')

xlswrite('Bio_sampling12.xlsx',fluxDist_rxns_R_12,'Sheet1','A1')
xlswrite('Bio_sampling12.xlsx',fluxDist_conditions_R_12,'Sheet1','B1')
xlswrite('Bio_sampling12.xlsx',fluxDist_R_12,'Sheet1','C1')

xlswrite('Bio_sampling13.xlsx',fluxDist_rxns_R_13,'Sheet1','A1')
xlswrite('Bio_sampling13.xlsx',fluxDist_conditions_R_13,'Sheet1','B1')
xlswrite('Bio_sampling13.xlsx',fluxDist_R_13,'Sheet1','C1')

xlswrite('Bio_sampling14.xlsx',fluxDist_rxns_R_14,'Sheet1','A1')
xlswrite('Bio_sampling14.xlsx',fluxDist_conditions_R_14,'Sheet1','B1')
xlswrite('Bio_sampling14.xlsx',fluxDist_R_14,'Sheet1','C1')

xlswrite('Bio_sampling15.xlsx',fluxDist_rxns_R_15,'Sheet1','A1')
xlswrite('Bio_sampling15.xlsx',fluxDist_conditions_R_15,'Sheet1','B1')
xlswrite('Bio_sampling15.xlsx',fluxDist_R_15,'Sheet1','C1')

xlswrite('Bio_sampling16.xlsx',fluxDist_rxns_R_16,'Sheet1','A1')
xlswrite('Bio_sampling16.xlsx',fluxDist_conditions_R_16,'Sheet1','B1')
xlswrite('Bio_sampling16.xlsx',fluxDist_R_16,'Sheet1','C1')

xlswrite('Bio_sampling17.xlsx',fluxDist_rxns_R_17,'Sheet1','A1')
xlswrite('Bio_sampling17.xlsx',fluxDist_conditions_R_17,'Sheet1','B1')
xlswrite('Bio_sampling17.xlsx',fluxDist_R_17,'Sheet1','C1')

xlswrite('Bio_sampling18.xlsx',fluxDist_rxns_R_18,'Sheet1','A1')
xlswrite('Bio_sampling18.xlsx',fluxDist_conditions_R_18,'Sheet1','B1')
xlswrite('Bio_sampling18.xlsx',fluxDist_R_18,'Sheet1','C1')

%% Reaction essentiality analysis on basal biomass models

[pa14Bio_essRxns, U5Bio_essRxns, U24Bio_essRxns, T5Bio_essRxns, T24Bio_essRxns] = persisterEssentialityRxnMax(pa14Bio, U5Bio, U24Bio, T5Bio, T24Bio);

U5Bio_essRxns = U5Bio_essRxns';
U24Bio_essRxns = U24Bio_essRxns';
T5Bio_essRxns = T5Bio_essRxns';
T24Bio_essRxns = T24Bio_essRxns';

%% Essentiality sensitivity analysis on basal biomass models
%[model_ess_pa14Bio, model_ess_U5Bio, model_ess_U24Bio, model_ess_T5Bio, model_ess_T24Bio, common_essBio, unique_treatedBio, unique_untreatedBio, unique_U5Bio, unique_U24Bio, unique_T5Bio, unique_T24Bio] = persisterEssentialityMin(pa14Bio, U5Bio, U24Bio, T5Bio, T24Bio);

pval_range = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
log2FC_range = [0:0.5:2.5];

flux_U5 = zeros(length(pval_range),length(log2FC_range));
flux_T5 = zeros(length(pval_range),length(log2FC_range));
flux_U24 = zeros(length(pval_range),length(log2FC_range));
flux_T24 = zeros(length(pval_range),length(log2FC_range));

metabsInt_U5 = zeros(length(pval_range),length(log2FC_range));
metabsInt_T5 = zeros(length(pval_range),length(log2FC_range));
metabsInt_U24 = zeros(length(pval_range),length(log2FC_range));
metabsInt_T24 = zeros(length(pval_range),length(log2FC_range));

essLength_U5 = zeros(length(pval_range),length(log2FC_range));
essLength_T5 = zeros(length(pval_range),length(log2FC_range));
essLength_U24 = zeros(length(pval_range),length(log2FC_range));
essLength_T24 = zeros(length(pval_range),length(log2FC_range));
essLength_common = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueTreated = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueUntreated = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueU5 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueT5 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueU24 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueT24 = zeros(length(pval_range),length(log2FC_range));

geneStatus_U5 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_T5 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_U24 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_T24 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));

geneStatusTotal_U5 = zeros(1,length(pa14.genes));
geneStatusTotal_T5 = zeros(1,length(pa14.genes));
geneStatusTotal_U24 = zeros(1,length(pa14.genes));
geneStatusTotal_T24 = zeros(1,length(pa14.genes));

for i = 1:length(pval_range)
    for j = 1:length(log2FC_range)
        
        pa14Bio = changeObjective(pa14,'PA14_Biomass',1);
        pa14Bio_sol = optimizeCbModel(pa14Bio);
        
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration_newmodel(pa14Bio, diffEX_file, RNA_pval, RNA_fc);

        U5_sol = optimizeCbModel(U5);
        T5_sol = optimizeCbModel(T5);
        U24_sol = optimizeCbModel(U24);
        T24_sol = optimizeCbModel(T24);

        flux_U5(i,j) = U5_sol.f;
        flux_T5(i,j) = T5_sol.f;
        flux_U24(i,j) = U24_sol.f;
        flux_T24(i,j) = T24_sol.f;
    
        metabsInt_U5(i,j) = length(untreated_metabs_5_int);
        metabsInt_T5(i,j) = length(persister_metabs_5_int);
        metabsInt_U24(i,j) = length(untreated_metabs_24_int);
        metabsInt_T24(i,j) = length(persister_metabs_24_int);
        
        U5.ub(1487) = 0.1*U5_sol.f;
        T5.ub(1487) = 0.1*T5_sol.f;
        U24.ub(1487) = 0.1*U24_sol.f;
        T24.ub(1487) = 0.1*T24_sol.f;
        pa14Bio.ub(1487) = 0.1*pa14Bio_sol.f;

        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMax(pa14, U5, U24, T5, T24);

        essLength_U5(i,j) = length(model_ess_U5);
        essLength_T5(i,j) = length(model_ess_T5);
        essLength_U24(i,j) = length(model_ess_U24);
        essLength_T24(i,j) = length(model_ess_T24);
        essLength_common(i,j) = length(common_ess);
        essLength_uniqueTreated(i,j) = length(unique_treated);
        essLength_uniqueUntreated(i,j) = length(unique_untreated);
        essLength_uniqueU5(i,j) = length(unique_U5);
        essLength_uniqueT5(i,j) = length(unique_T5);
        essLength_uniqueU24(i,j) = length(unique_U24);
        essLength_uniqueT24(i,j) = length(unique_T24);
    
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomas_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

        for k = 1:length(pa14.genes)
            geneTest_U5 = isempty(find(not(cellfun('isempty',strfind(model_ess_U5, pa14.genes(k))))));
            geneTest_T5 = isempty(find(not(cellfun('isempty',strfind(model_ess_T5, pa14.genes(k))))));
            geneTest_U24 = isempty(find(not(cellfun('isempty',strfind(model_ess_U24, pa14.genes(k))))));
            geneTest_T24 = isempty(find(not(cellfun('isempty',strfind(model_ess_T24, pa14.genes(k))))));

            if geneTest_U5 == 0
                geneStatus_U5(i,j,k) = geneStatus_U5(i,j,k) + 1;
                geneStatusTotal_U5(k) = geneStatusTotal_U5(k) + 1;
            end
            if geneTest_T5 == 0
                geneStatus_T5(i,j,k) = geneStatus_T5(i,j,k) + 1;
                geneStatusTotal_T5(k) = geneStatusTotal_T5(k) + 1;
            end
            if geneTest_U24 == 0
                geneStatus_U24(i,j,k) = geneStatus_U24(i,j,k) + 1;
                geneStatusTotal_U24(k) = geneStatusTotal_U24(k) + 1;
            end
            if geneTest_T24 == 0
                geneStatus_T24(i,j,k) = geneStatus_T24(i,j,k) + 1;
                geneStatusTotal_T24(k) = geneStatusTotal_T24(k) + 1;
            end
        end 
    end
end

xlswrite('180921_persisterSensitivityAnalysis_Biomass_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
xlswrite('180921_persisterSensitivityAnalysis_Biomass_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
xlswrite('180921_persisterSensitivityAnalysis_Biomass_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
xlswrite('180921_persisterSensitivityAnalysis_Biomass_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
xlswrite('180921_persisterSensitivityAnalysis_Biomass_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');












%% Create basal ATP models
pa14ATP = addDemandReaction(pa14, 'cpd00002[cytosol]');
pa14ATP = changeObjective(pa14ATP, 'DM_cpd00002[cytosol]',1);
pa14ATP_sol = optimizeCbModel(pa14ATP);
pa14ATP = creategrRulesField(pa14ATP);

RNA_pval = [0.5];
RNA_fc = [0];

[U5ATP, U24ATP, T5ATP, T24ATP, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration_newmodel(pa14ATP, diffEX_file, RNA_pval, RNA_fc);

U5ATP_sol = optimizeCbModel(U5ATP);
T5ATP_sol = optimizeCbModel(T5ATP);
U24ATP_sol = optimizeCbModel(U24ATP);
T24ATP_sol = optimizeCbModel(T24ATP);

% need to modify this
U5ATP.ub(length(pa14ATP.rxns)) = 0.1*U5ATP_sol.f;
T5ATP.ub(length(pa14ATP.rxns)) = 0.1*T5ATP_sol.f;
U24ATP.ub(length(pa14ATP.rxns)) = 0.1*U24ATP_sol.f;
T24ATP.ub(length(pa14ATP.rxns)) = 0.1*T24ATP_sol.f;
pa14ATP.ub(length(pa14ATP.rxns)) = 0.1*pa14ATP_sol.f;

%% Flux sampling on basal ATP models
% Note, in order to run this next chunk, I needed to change the directory 
% to optGpSampler_1.1_Matlab

flux_U5ATP = optGpSampler(U5ATP,[],3000,1,1,'gurobi',0);
flux_U24ATP = optGpSampler(U24ATP,[],3000,1,1,'gurobi',0);
flux_T5ATP = optGpSampler(T5ATP,[],3000,1,1,'gurobi',0);
flux_T24ATP = optGpSampler(T24ATP,[],3000,1,1,'gurobi',0);

fluxDist_U5ATP = zeros(length(U5ATP.rxns),3000);
fluxMean_U5ATP = zeros(length(U5ATP.rxns),1);
fluxMedian_U5ATP = zeros(length(U5ATP.rxns),1);
for i = 1:length(U5ATP.rxns)
    fluxDist_U5ATP(i,:) = flux_U5ATP.points(findRxnIDs(U5ATP,U5ATP.rxns(i)),:);
    fluxMean_U5ATP(i) = mean(fluxDist_U5ATP(i,:));
    fluxMedian_U5ATP(i) = median(fluxDist_U5ATP(i,:));
end

fluxDist_U24ATP = zeros(length(U24ATP.rxns),3000);
fluxMean_U24ATP = zeros(length(U24ATP.rxns),1);
fluxMedian_U24ATP = zeros(length(U24ATP.rxns),1);
for i = 1:length(U24ATP.rxns)
    fluxDist_U24ATP(i,:) = flux_U24ATP.points(findRxnIDs(U24ATP,U24ATP.rxns(i)),:);
    fluxMean_U24ATP(i) = mean(fluxDist_U24ATP(i,:));
    fluxMedian_U24ATP(i) = median(fluxDist_U24ATP(i,:));
end

fluxDist_T5ATP = zeros(length(T5ATP.rxns),3000);
fluxMean_T5ATP = zeros(length(T5ATP.rxns),1);
fluxMedian_T5ATP = zeros(length(T5ATP.rxns),1);
for i = 1:length(T5ATP.rxns)
    fluxDist_T5ATP(i,:) = flux_T5ATP.points(findRxnIDs(T5ATP,T5ATP.rxns(i)),:);
    fluxMean_T5ATP(i) = mean(fluxDist_T5ATP(i,:));
    fluxMedian_T5ATP(i) = median(fluxDist_T5ATP(i,:));
end

fluxDist_T24ATP = zeros(length(T24ATP.rxns),3000);
fluxMean_T24ATP = zeros(length(T24ATP.rxns),1);
fluxMedian_T24ATP = zeros(length(T24ATP.rxns),1);
for i = 1:length(T24ATP.rxns)
    fluxDist_T24ATP(i,:) = flux_T24ATP.points(findRxnIDs(T24ATP,T24ATP.rxns(i)),:);
    fluxMean_T24ATP(i) = mean(fluxDist_T24ATP(i,:));
    fluxMedian_T24ATP(i) = median(fluxDist_T24ATP(i,:));
end

% for each reaction, find the proportion of fluxes greater than 0.00001
proportion = zeros(length(pa14.rxns),4);

for i = 1:length(pa14.rxns)
    activeSamples_U5ATP = find(fluxDist_U5ATP(i,:) > 0.0001);
    activeSamples_U24ATP = find(fluxDist_U24ATP(i,:) > 0.0001);
    activeSamples_T5ATP = find(fluxDist_T5ATP(i,:) > 0.0001);
    activeSamples_T24ATP = find(fluxDist_T24ATP(i,:) > 0.0001);

    proportion(i,1) = length(activeSamples_U5ATP)/length(fluxDist_U5ATP);
    proportion(i,2) = length(activeSamples_U24ATP)/length(fluxDist_U24ATP);
    proportion(i,3) = length(activeSamples_T5ATP)/length(fluxDist_T5ATP);
    proportion(i,4) = length(activeSamples_T24ATP)/length(fluxDist_T24ATP);
end

% tidy data for R

fluxDist_U5ATP_R = [];
fluxDist_U24ATP_R = [];
fluxDist_T5ATP_R = [];
fluxDist_T24ATP_R = [];
fluxDistATP = [];

for i = 1:length(pa14.rxns)
    fluxDist_U5ATP_R = [fluxDist_U5ATP(i,:)'];
    fluxDist_U24ATP_R = [fluxDist_U24ATP(i,:)'];
    fluxDist_T5ATP_R = [fluxDist_T5ATP(i,:)'];
    fluxDist_T24ATP_R = [fluxDist_T24ATP(i,:)'];
    fluxDistATP = [fluxDistATP;
                fluxDist_U5ATP_R;
                fluxDist_U24ATP_R;
                fluxDist_T5ATP_R;
                fluxDist_T24ATP_R];
end

fluxDist_rxns_R = [];
fluxDist_conditions_R = [];

conditions_R = [repmat({'U5'},[3000 1]);
                repmat({'U24'},[3000 1]);
                repmat({'T5'}, [3000 1]);
                repmat({'T24'}, [3000 1])];

for i = 1:length(pa14.rxns)
    rxns_mat = repmat(pa14.rxns(i),(3000*4),1);
    if i == 1
        fluxDist_rxns_R = rxns_mat;
    else
        fluxDist_rxns_R = [fluxDist_rxns_R;rxns_mat];
    end
    
    fluxDist_conditions_R = [fluxDist_conditions_R; conditions_R];
end

% I can't save these vectors in excel cause they are too large.
% The max length is 1048576

parse = 17940000/1048576;

fluxDist_rxns_R_1 = fluxDist_rxns_R(1:1048576);
fluxDist_rxns_R_2 = fluxDist_rxns_R(1048577:2097152);
fluxDist_rxns_R_3 = fluxDist_rxns_R(2097153:3145728);
fluxDist_rxns_R_4 = fluxDist_rxns_R(3145729:4194300);
fluxDist_rxns_R_5 = fluxDist_rxns_R(4194301:5242871);
fluxDist_rxns_R_6 = fluxDist_rxns_R(5242872:6291447);
fluxDist_rxns_R_7 = fluxDist_rxns_R(6291448:7340023);
fluxDist_rxns_R_8 = fluxDist_rxns_R(7340024:8388509);
fluxDist_rxns_R_9 = fluxDist_rxns_R(8388510:9437085);
fluxDist_rxns_R_10 = fluxDist_rxns_R(9437086:10485661);
fluxDist_rxns_R_11 = fluxDist_rxns_R(10485662:11534237);
fluxDist_rxns_R_12 = fluxDist_rxns_R(11534238:12582807);
fluxDist_rxns_R_13 = fluxDist_rxns_R(12582808:13631383);
fluxDist_rxns_R_14 = fluxDist_rxns_R(13631384:14679959);
fluxDist_rxns_R_15 = fluxDist_rxns_R(14679960:15728535);
fluxDist_rxns_R_16 = fluxDist_rxns_R(15728536:16777111);
fluxDist_rxns_R_17 = fluxDist_rxns_R(16777112:17825687);
fluxDist_rxns_R_18 = fluxDist_rxns_R(17825688:end);

fluxDist_conditions_R_1 = fluxDist_conditions_R(1:1048576);
fluxDist_conditions_R_2 = fluxDist_conditions_R(1048577:2097152);
fluxDist_conditions_R_3 = fluxDist_conditions_R(2097153:3145728);
fluxDist_conditions_R_4 = fluxDist_conditions_R(3145729:4194300);
fluxDist_conditions_R_5 = fluxDist_conditions_R(4194301:5242871);
fluxDist_conditions_R_6 = fluxDist_conditions_R(5242872:6291447);
fluxDist_conditions_R_7 = fluxDist_conditions_R(6291448:7340023);
fluxDist_conditions_R_8 = fluxDist_conditions_R(7340024:8388509);
fluxDist_conditions_R_9 = fluxDist_conditions_R(8388510:9437085);
fluxDist_conditions_R_10 = fluxDist_conditions_R(9437086:10485661);
fluxDist_conditions_R_11 = fluxDist_conditions_R(10485662:11534237);
fluxDist_conditions_R_12 = fluxDist_conditions_R(11534238:12582807);
fluxDist_conditions_R_13 = fluxDist_conditions_R(12582808:13631383);
fluxDist_conditions_R_14 = fluxDist_conditions_R(13631384:14679959);
fluxDist_conditions_R_15 = fluxDist_conditions_R(14679960:15728535);
fluxDist_conditions_R_16 = fluxDist_conditions_R(15728536:16777111);
fluxDist_conditions_R_17 = fluxDist_conditions_R(16777112:17825687);
fluxDist_conditions_R_18 = fluxDist_conditions_R(17825688:end);

fluxDist_R_1 = fluxDistATP(1:1048576);
fluxDist_R_2 = fluxDistATP(1048577:2097152);
fluxDist_R_3 = fluxDistATP(2097153:3145728);
fluxDist_R_4 = fluxDistATP(3145729:4194300);
fluxDist_R_5 = fluxDistATP(4194301:5242871);
fluxDist_R_6 = fluxDistATP(5242872:6291447);
fluxDist_R_7 = fluxDistATP(6291448:7340023);
fluxDist_R_8 = fluxDistATP(7340024:8388509);
fluxDist_R_9 = fluxDistATP(8388510:9437085);
fluxDist_R_10 = fluxDistATP(9437086:10485661);
fluxDist_R_11 = fluxDistATP(10485662:11534237);
fluxDist_R_12 = fluxDistATP(11534238:12582807);
fluxDist_R_13 = fluxDistATP(12582808:13631383);
fluxDist_R_14 = fluxDistATP(13631384:14679959);
fluxDist_R_15 = fluxDistATP(14679960:15728535);
fluxDist_R_16 = fluxDistATP(15728536:16777111);
fluxDist_R_17 = fluxDistATP(16777112:17825687);
fluxDist_R_18 = fluxDistATP(17825688:end);

xlswrite('ATP_sampling1.xlsx',fluxDist_rxns_R_1,'Sheet1','A1')
xlswrite('ATP_sampling1.xlsx',fluxDist_conditions_R_1,'Sheet1','B1')
xlswrite('ATP_sampling1.xlsx',fluxDist_R_1,'Sheet1','C1')

xlswrite('ATP_sampling2.xlsx',fluxDist_rxns_R_2,'Sheet1','A1')
xlswrite('ATP_sampling2.xlsx',fluxDist_conditions_R_2,'Sheet1','B1')
xlswrite('ATP_sampling2.xlsx',fluxDist_R_2,'Sheet1','C1')

xlswrite('ATP_sampling3.xlsx',fluxDist_rxns_R_3,'Sheet1','A1')
xlswrite('ATP_sampling3.xlsx',fluxDist_conditions_R_3,'Sheet1','B1')
xlswrite('ATP_sampling3.xlsx',fluxDist_R_3,'Sheet1','C1')

xlswrite('ATP_sampling4.xlsx',fluxDist_rxns_R_4,'Sheet1','A1')
xlswrite('ATP_sampling4.xlsx',fluxDist_conditions_R_4,'Sheet1','B1')
xlswrite('ATP_sampling4.xlsx',fluxDist_R_4,'Sheet1','C1')

xlswrite('ATP_sampling5.xlsx',fluxDist_rxns_R_5,'Sheet1','A1')
xlswrite('ATP_sampling5.xlsx',fluxDist_conditions_R_5,'Sheet1','B1')
xlswrite('ATP_sampling5.xlsx',fluxDist_R_5,'Sheet1','C1')

xlswrite('ATP_sampling6.xlsx',fluxDist_rxns_R_6,'Sheet1','A1')
xlswrite('ATP_sampling6.xlsx',fluxDist_conditions_R_6,'Sheet1','B1')
xlswrite('ATP_sampling6.xlsx',fluxDist_R_6,'Sheet1','C1')

xlswrite('ATP_sampling7.xlsx',fluxDist_rxns_R_7,'Sheet1','A1')
xlswrite('ATP_sampling7.xlsx',fluxDist_conditions_R_7,'Sheet1','B1')
xlswrite('ATP_sampling7.xlsx',fluxDist_R_7,'Sheet1','C1')

xlswrite('ATP_sampling8.xlsx',fluxDist_rxns_R_8,'Sheet1','A1')
xlswrite('ATP_sampling8.xlsx',fluxDist_conditions_R_8,'Sheet1','B1')
xlswrite('ATP_sampling8.xlsx',fluxDist_R_8,'Sheet1','C1')

xlswrite('ATP_sampling9.xlsx',fluxDist_rxns_R_9,'Sheet1','A1')
xlswrite('ATP_sampling9.xlsx',fluxDist_conditions_R_9,'Sheet1','B1')
xlswrite('ATP_sampling9.xlsx',fluxDist_R_9,'Sheet1','C1')

xlswrite('ATP_sampling10.xlsx',fluxDist_rxns_R_10,'Sheet1','A1')
xlswrite('ATP_sampling10.xlsx',fluxDist_conditions_R_10,'Sheet1','B1')
xlswrite('ATP_sampling10.xlsx',fluxDist_R_10,'Sheet1','C1')

xlswrite('ATP_sampling11.xlsx',fluxDist_rxns_R_11,'Sheet1','A1')
xlswrite('ATP_sampling11.xlsx',fluxDist_conditions_R_11,'Sheet1','B1')
xlswrite('ATP_sampling11.xlsx',fluxDist_R_11,'Sheet1','C1')

xlswrite('ATP_sampling12.xlsx',fluxDist_rxns_R_12,'Sheet1','A1')
xlswrite('ATP_sampling12.xlsx',fluxDist_conditions_R_12,'Sheet1','B1')
xlswrite('ATP_sampling12.xlsx',fluxDist_R_12,'Sheet1','C1')

xlswrite('ATP_sampling13.xlsx',fluxDist_rxns_R_13,'Sheet1','A1')
xlswrite('ATP_sampling13.xlsx',fluxDist_conditions_R_13,'Sheet1','B1')
xlswrite('ATP_sampling13.xlsx',fluxDist_R_13,'Sheet1','C1')

xlswrite('ATP_sampling14.xlsx',fluxDist_rxns_R_14,'Sheet1','A1')
xlswrite('ATP_sampling14.xlsx',fluxDist_conditions_R_14,'Sheet1','B1')
xlswrite('ATP_sampling14.xlsx',fluxDist_R_14,'Sheet1','C1')

xlswrite('ATP_sampling15.xlsx',fluxDist_rxns_R_15,'Sheet1','A1')
xlswrite('ATP_sampling15.xlsx',fluxDist_conditions_R_15,'Sheet1','B1')
xlswrite('ATP_sampling15.xlsx',fluxDist_R_15,'Sheet1','C1')

xlswrite('ATP_sampling16.xlsx',fluxDist_rxns_R_16,'Sheet1','A1')
xlswrite('ATP_sampling16.xlsx',fluxDist_conditions_R_16,'Sheet1','B1')
xlswrite('ATP_sampling16.xlsx',fluxDist_R_16,'Sheet1','C1')

xlswrite('ATP_sampling17.xlsx',fluxDist_rxns_R_17,'Sheet1','A1')
xlswrite('ATP_sampling17.xlsx',fluxDist_conditions_R_17,'Sheet1','B1')
xlswrite('ATP_sampling17.xlsx',fluxDist_R_17,'Sheet1','C1')

xlswrite('ATP_sampling18.xlsx',fluxDist_rxns_R_18,'Sheet1','A1')
xlswrite('ATP_sampling18.xlsx',fluxDist_conditions_R_18,'Sheet1','B1')
xlswrite('ATP_sampling18.xlsx',fluxDist_R_18,'Sheet1','C1')

%% Reaction essentiality analysis on basal ATP models

[pa14ATP_essRxns, U5ATP_essRxns, U24ATP_essRxns, T5ATP_essRxns, T24ATP_essRxns] = persisterEssentialityRxnMax(pa14ATP, U5ATP, U24ATP, T5ATP, T24ATP);

U5ATP_essRxns = U5ATP_essRxns';
T5ATP_essRxns = T5ATP_essRxns';
U24ATP_essRxns = U24ATP_essRxns';
T24ATP_essRxns = T24ATP_essRxns';

%% Essentiality analysis on basal ATP models

pval_range = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
log2FC_range = [0:0.5:2.5];

flux_U5 = zeros(length(pval_range),length(log2FC_range));
flux_T5 = zeros(length(pval_range),length(log2FC_range));
flux_U24 = zeros(length(pval_range),length(log2FC_range));
flux_T24 = zeros(length(pval_range),length(log2FC_range));

metabsInt_U5 = zeros(length(pval_range),length(log2FC_range));
metabsInt_T5 = zeros(length(pval_range),length(log2FC_range));
metabsInt_U24 = zeros(length(pval_range),length(log2FC_range));
metabsInt_T24 = zeros(length(pval_range),length(log2FC_range));

essLength_U5 = zeros(length(pval_range),length(log2FC_range));
essLength_T5 = zeros(length(pval_range),length(log2FC_range));
essLength_U24 = zeros(length(pval_range),length(log2FC_range));
essLength_T24 = zeros(length(pval_range),length(log2FC_range));
essLength_common = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueTreated = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueUntreated = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueU5 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueT5 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueU24 = zeros(length(pval_range),length(log2FC_range));
essLength_uniqueT24 = zeros(length(pval_range),length(log2FC_range));

geneStatus_U5 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_T5 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_U24 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));
geneStatus_T24 = zeros(length(pval_range),length(log2FC_range),length(pa14.genes));

geneStatusTotal_U5 = zeros(1,length(pa14.genes));
geneStatusTotal_T5 = zeros(1,length(pa14.genes));
geneStatusTotal_U24 = zeros(1,length(pa14.genes));
geneStatusTotal_T24 = zeros(1,length(pa14.genes));

for i = 1:length(pval_range)
    for j = 1:length(log2FC_range)
        
        pa14ATP = addDemandReaction(pa14, 'cpd00002[c]');
        pa14ATP = changeObjective(pa14ATP, 'DM_cpd00002[c]',1);
        pa14ATP_sol = optimizeCbModel(pa14ATP);
        
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration_newmodel(pa14ATP, diffEX_file, RNA_pval, RNA_fc);

        U5_sol = optimizeCbModel(U5);
        T5_sol = optimizeCbModel(T5);
        U24_sol = optimizeCbModel(U24);
        T24_sol = optimizeCbModel(T24);

        flux_U5(i,j) = U5_sol.f;
        flux_T5(i,j) = T5_sol.f;
        flux_U24(i,j) = U24_sol.f;
        flux_T24(i,j) = T24_sol.f;
    
        metabsInt_U5(i,j) = length(untreated_metabs_5_int);
        metabsInt_T5(i,j) = length(persister_metabs_5_int);
        metabsInt_U24(i,j) = length(untreated_metabs_24_int);
        metabsInt_T24(i,j) = length(persister_metabs_24_int);
        
        U5.ub(length(pa14ATP.rxns)) = 0.1*U5_sol.f;
        T5.ub(length(pa14ATP.rxns)) = 0.1*T5_sol.f;
        U24.ub(length(pa14ATP.rxns)) = 0.1*U24_sol.f;
        T24.ub(length(pa14ATP.rxns)) = 0.1*T24_sol.f;
        pa14ATP.ub(length(pa14ATP.rxns)) = 0.1*pa14ATP_sol.f;

        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMax(pa14, U5, U24, T5, T24);

        essLength_U5(i,j) = length(model_ess_U5);
        essLength_T5(i,j) = length(model_ess_T5);
        essLength_U24(i,j) = length(model_ess_U24);
        essLength_T24(i,j) = length(model_ess_T24);
        essLength_common(i,j) = length(common_ess);
        essLength_uniqueTreated(i,j) = length(unique_treated);
        essLength_uniqueUntreated(i,j) = length(unique_untreated);
        essLength_uniqueU5(i,j) = length(unique_U5);
        essLength_uniqueT5(i,j) = length(unique_T5);
        essLength_uniqueU24(i,j) = length(unique_U24);
        essLength_uniqueT24(i,j) = length(unique_T24);
    
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinBiomas_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

        for k = 1:length(pa14.genes)
            geneTest_U5 = isempty(find(not(cellfun('isempty',strfind(model_ess_U5, pa14.genes(k))))));
            geneTest_T5 = isempty(find(not(cellfun('isempty',strfind(model_ess_T5, pa14.genes(k))))));
            geneTest_U24 = isempty(find(not(cellfun('isempty',strfind(model_ess_U24, pa14.genes(k))))));
            geneTest_T24 = isempty(find(not(cellfun('isempty',strfind(model_ess_T24, pa14.genes(k))))));

            if geneTest_U5 == 0
                geneStatus_U5(i,j,k) = geneStatus_U5(i,j,k) + 1;
                geneStatusTotal_U5(k) = geneStatusTotal_U5(k) + 1;
            end
            if geneTest_T5 == 0
                geneStatus_T5(i,j,k) = geneStatus_T5(i,j,k) + 1;
                geneStatusTotal_T5(k) = geneStatusTotal_T5(k) + 1;
            end
            if geneTest_U24 == 0
                geneStatus_U24(i,j,k) = geneStatus_U24(i,j,k) + 1;
                geneStatusTotal_U24(k) = geneStatusTotal_U24(k) + 1;
            end
            if geneTest_T24 == 0
                geneStatus_T24(i,j,k) = geneStatus_T24(i,j,k) + 1;
                geneStatusTotal_T24(k) = geneStatusTotal_T24(k) + 1;
            end
        end 
    end
end

xlswrite('180921_persisterSensitivityAnalysis_ATP_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
xlswrite('180921_persisterSensitivityAnalysis_ATP_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
xlswrite('180921_persisterSensitivityAnalysis_ATP_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
xlswrite('180921_persisterSensitivityAnalysis_ATP_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
xlswrite('180921_persisterSensitivityAnalysis_ATP_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');
