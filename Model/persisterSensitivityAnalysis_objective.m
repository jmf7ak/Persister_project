%% initialize environment

clear all
close all
clc
initCobraToolbox()
changeCobraSolver('gurobi');

% Windows:
 pa14 = xls2model_JAB('C:\Users\jmfic\OneDrive\Documents\21_SUMMA\Persister\data\PA14recon1_v24_published.xlsx');
% Mac:
% pa14 = d_xls2model_JAB('/Users/papinlab/Dropbox/pseudomonas_model/PA14recon1_v24_published.xlsx');

pa14 = changeMedia_SEED(pa14, 1,'');
pa14 = changeObjective(pa14,'PA14_Biomass',1);
pa14_sol = optimizeCbModel(pa14);

% Windows:
diffEX_file = 'C:\Users\jmfic\OneDrive\Documents\21_SUMMA\Persister\data\DESeq_dds_groupbatch_all2metab_integration_180221.xls';
% Mac:
% diffEX_file = '/Users/papinlab/Dropbox/Lab/Projects/Persisters/RNA_Sequencing/BIT_RNA_Sequencing/Second-run/Analysis/DESeq_dds_groupbatch_all2metab_integration_180221.xls';


%% Pval Sensitivity Analysis

pval_range = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];

flux_U5 = zeros(1,length(pval_range));
flux_T5 = zeros(1,length(pval_range));
flux_U24 = zeros(1,length(pval_range));
flux_T24 = zeros(1,length(pval_range));

metabsInt_U5 = zeros(1,length(pval_range));
metabsInt_T5 = zeros(1,length(pval_range));
metabsInt_U24 = zeros(1,length(pval_range));
metabsInt_T24 = zeros(1,length(pval_range));

essLength_U5 = zeros(1,length(pval_range));
essLength_T5 = zeros(1,length(pval_range));
essLength_U24 = zeros(1,length(pval_range));
essLength_T24 = zeros(1,length(pval_range));
essLength_common = zeros(1,length(pval_range));
essLength_uniqueTreated = zeros(1,length(pval_range));
essLength_uniqueUntreated = zeros(1,length(pval_range));
essLength_uniqueU5 = zeros(1,length(pval_range));
essLength_uniqueT5 = zeros(1,length(pval_range));
essLength_uniqueU24 = zeros(1,length(pval_range));
essLength_uniqueT24 = zeros(1,length(pval_range));

geneStatus_U5 = zeros(1,length(pa14.genes));
geneStatus_T5 = zeros(1,length(pa14.genes));
geneStatus_U24 = zeros(1,length(pa14.genes));
geneStatus_T24 = zeros(1,length(pa14.genes));

for i = 1:length(pval_range)
    RNA_pval = pval_range(i);
    RNA_fc = 0;
    
    [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);

    U5_sol = optimizeCbModel(U5);
    T5_sol = optimizeCbModel(T5);
    U24_sol = optimizeCbModel(U24);
    T24_sol = optimizeCbModel(T24);

    flux_U5(i) = U5_sol.f;
    flux_T5(i) = T5_sol.f;
    flux_U24(i) = U24_sol.f;
    flux_T24(i) = T24_sol.f;
    
    metabsInt_U5(i) = length(untreated_metabs_5_int);
    metabsInt_T5(i) = length(persister_metabs_5_int);
    metabsInt_U24(i) = length(untreated_metabs_24_int);
    metabsInt_T24(i) = length(persister_metabs_24_int);
    
    U5.lb(1487) = 0.1*U5_sol.f;
    T5.lb(1487) = 0.1*T5_sol.f;
    U24.lb(1487) = 0.1*U24_sol.f;
    T24.lb(1487) = 0.1*T24_sol.f;
    pa14.lb(1487) = 0.1*pa14_sol.f;

    [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMin(pa14, U5, U24, T5, T24);

    essLength_U5(i) = length(model_ess_U5);
    essLength_T5(i) = length(model_ess_T5);
    essLength_U24(i) = length(model_ess_U24);
    essLength_T24(i) = length(model_ess_T24);
    essLength_common(i) = length(common_ess);
    essLength_uniqueTreated(i) = length(unique_treated);
    essLength_uniqueUntreated(i) = length(unique_untreated);
    essLength_uniqueU5(i) = length(unique_U5);
    essLength_uniqueT5(i) = length(unique_T5);
    essLength_uniqueU24(i) = length(unique_U24);
    essLength_uniqueT24(i) = length(unique_T24);
    
%     xlswrite('180723_persisterSensitivityAnalysis_pval_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval)),'A1');
%     xlswrite('180723_persisterSensitivityAnalysis_pval_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval)),'B1');
%     xlswrite('180723_persisterSensitivityAnalysis_pval_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval)),'C1');
%     xlswrite('180723_persisterSensitivityAnalysis_pval_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval)),'D1');

    for j = 1:length(pa14.genes)
        geneTest_U5 = isempty(find(not(cellfun('isempty',strfind(model_ess_U5, pa14.genes(j))))));
        geneTest_T5 = isempty(find(not(cellfun('isempty',strfind(model_ess_T5, pa14.genes(j))))));
        geneTest_U24 = isempty(find(not(cellfun('isempty',strfind(model_ess_U24, pa14.genes(j))))));
        geneTest_T24 = isempty(find(not(cellfun('isempty',strfind(model_ess_T24, pa14.genes(j))))));

        if geneTest_U5 == 0
            geneStatus_U5(j) = geneStatus_U5(j) + 1;
        end
        if geneTest_T5 == 0
            geneStatus_T5(j) = geneStatus_T5(j) + 1;
        end
        if geneTest_U24 == 0
            geneStatus_U24(j) = geneStatus_U24(j) + 1;
        end
        if geneTest_T24 == 0
            geneStatus_T24(j) = geneStatus_T24(j) + 1;
        end
    end 
end

% xlswrite('180723_persisterSensitivityAnalysis_pval_geneStatus.xlsx',pa14.genes,'Sheet1','A1');
% xlswrite('180723_persisterSensitivityAnalysis_pval_geneStatus.xlsx',geneStatus_U5','Sheet1','B1');
% xlswrite('180723_persisterSensitivityAnalysis_pval_geneStatus.xlsx',geneStatus_T5','Sheet1','C1');
% xlswrite('180723_persisterSensitivityAnalysis_pval_geneStatus.xlsx',geneStatus_U24','Sheet1','D1');
% xlswrite('180723_persisterSensitivityAnalysis_pval_geneStatus.xlsx',geneStatus_T24','Sheet1','E1');

%% log2FC Sensitivity Analysis

log2FC_range = [0:0.5:2.5];

flux_U5 = zeros(1,length(log2FC_range));
flux_T5 = zeros(1,length(log2FC_range));
flux_U24 = zeros(1,length(log2FC_range));
flux_T24 = zeros(1,length(log2FC_range));

metabsInt_U5 = zeros(1,length(log2FC_range));
metabsInt_T5 = zeros(1,length(log2FC_range));
metabsInt_U24 = zeros(1,length(log2FC_range));
metabsInt_T24 = zeros(1,length(log2FC_range));

essLength_U5 = zeros(1,length(log2FC_range));
essLength_T5 = zeros(1,length(log2FC_range));
essLength_U24 = zeros(1,length(log2FC_range));
essLength_T24 = zeros(1,length(log2FC_range));
essLength_common = zeros(1,length(log2FC_range));
essLength_uniqueTreated = zeros(1,length(log2FC_range));
essLength_uniqueUntreated = zeros(1,length(log2FC_range));
essLength_uniqueU5 = zeros(1,length(log2FC_range));
essLength_uniqueT5 = zeros(1,length(log2FC_range));
essLength_uniqueU24 = zeros(1,length(log2FC_range));
essLength_uniqueT24 = zeros(1,length(log2FC_range));

geneStatus_U5 = zeros(1,length(pa14.genes));
geneStatus_T5 = zeros(1,length(pa14.genes));
geneStatus_U24 = zeros(1,length(pa14.genes));
geneStatus_T24 = zeros(1,length(pa14.genes));

for i = 1:length(log2FC_range)
    RNA_pval = 0.5;
    RNA_fc = log2FC_range(i);
    
    [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);

    U5_sol = optimizeCbModel(U5);
    T5_sol = optimizeCbModel(T5);
    U24_sol = optimizeCbModel(U24);
    T24_sol = optimizeCbModel(T24);

    flux_U5(i) = U5_sol.f;
    flux_T5(i) = T5_sol.f;
    flux_U24(i) = U24_sol.f;
    flux_T24(i) = T24_sol.f;
    
    metabsInt_U5(i) = length(untreated_metabs_5_int);
    metabsInt_T5(i) = length(persister_metabs_5_int);
    metabsInt_U24(i) = length(untreated_metabs_24_int);
    metabsInt_T24(i) = length(persister_metabs_24_int);

    [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentiality(pa14, U5, U24, T5, T24);

    essLength_U5(i) = length(model_ess_U5);
    essLength_T5(i) = length(model_ess_T5);
    essLength_U24(i) = length(model_ess_U24);
    essLength_T24(i) = length(model_ess_T24);
    essLength_common(i) = length(common_ess);
    essLength_uniqueTreated(i) = length(unique_treated);
    essLength_uniqueUntreated(i) = length(unique_untreated);
    essLength_uniqueU5(i) = length(unique_U5);
    essLength_uniqueT5(i) = length(unique_T5);
    essLength_uniqueU24(i) = length(unique_U24);
    essLength_uniqueT24(i) = length(unique_T24);
    
    xlswrite('180723_persisterSensitivityAnalysis_FC_essGenes.xlsx',model_ess_U5,strcat('log2FC=',num2str(RNA_fc)),'A1');
    xlswrite('180723_persisterSensitivityAnalysis_FC_essGenes.xlsx',model_ess_T5,strcat('log2FC=',num2str(RNA_fc)),'B1');
    xlswrite('180723_persisterSensitivityAnalysis_FC_essGenes.xlsx',model_ess_U24,strcat('log2FC=',num2str(RNA_fc)),'C1');
    xlswrite('180723_persisterSensitivityAnalysis_FC_essGenes.xlsx',model_ess_T24,strcat('log2FC=',num2str(RNA_fc)),'D1');

    for j = 1:length(pa14.genes)
        geneTest_U5 = isempty(find(not(cellfun('isempty',strfind(model_ess_U5, pa14.genes(j))))));
        geneTest_T5 = isempty(find(not(cellfun('isempty',strfind(model_ess_T5, pa14.genes(j))))));
        geneTest_U24 = isempty(find(not(cellfun('isempty',strfind(model_ess_U24, pa14.genes(j))))));
        geneTest_T24 = isempty(find(not(cellfun('isempty',strfind(model_ess_T24, pa14.genes(j))))));

        if geneTest_U5 == 0
            geneStatus_U5(j) = geneStatus_U5(j) + 1;
        end
        if geneTest_T5 == 0
            geneStatus_T5(j) = geneStatus_T5(j) + 1;
        end
        if geneTest_U24 == 0
            geneStatus_U24(j) = geneStatus_U24(j) + 1;
        end
        if geneTest_T24 == 0
            geneStatus_T24(j) = geneStatus_T24(j) + 1;
        end
    end 
    
end

xlswrite('180723_persisterSensitivityAnalysis_FC_geneStatus.xlsx',pa14.genes,'Sheet1','A1');
xlswrite('180723_persisterSensitivityAnalysis_FC_geneStatus.xlsx',geneStatus_U5','Sheet1','B1');
xlswrite('180723_persisterSensitivityAnalysis_FC_geneStatus.xlsx',geneStatus_T5','Sheet1','C1');
xlswrite('180723_persisterSensitivityAnalysis_FC_geneStatus.xlsx',geneStatus_U24','Sheet1','D1');
xlswrite('180723_persisterSensitivityAnalysis_FC_geneStatus.xlsx',geneStatus_T24','Sheet1','E1');

%% Sensitivity Analysis - minimize biomass objective

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
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);

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
        
        U5.lb(1487) = 0.1*U5_sol.f;
        T5.lb(1487) = 0.1*T5_sol.f;
        U24.lb(1487) = 0.1*U24_sol.f;
        T24.lb(1487) = 0.1*T24_sol.f;
        pa14.lb(1487) = 0.1*pa14_sol.f;

        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMin(pa14, U5, U24, T5, T24);

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
    
        xlswrite('181119_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
        xlswrite('181119_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
        xlswrite('181119_persisterSensitivityAnalysis_MinBiomass_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
        xlswrite('181119_persisterSensitivityAnalysis_MinBiomas_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

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

  xlswrite('180911_persisterSensitivityAnalysis_minBiomass_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
  xlswrite('180911_persisterSensitivityAnalysis_minBiomass_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
  xlswrite('180911_persisterSensitivityAnalysis_minBiomass_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
  xlswrite('180911_persisterSensitivityAnalysis_minBiomass_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
  xlswrite('180911_persisterSensitivityAnalysis_minBiomass_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');

  %% Sensitivity Analysis - minimize ATP production while forcing biomass

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
        
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);
        
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
        
        U5.lb(1487) = 0.1*U5_sol.f;
        T5.lb(1487) = 0.1*T5_sol.f;
        U24.lb(1487) = 0.1*U24_sol.f;
        T24.lb(1487) = 0.1*T24_sol.f;
        pa14.lb(1487) = 0.1*pa14_sol.f;
        
        U5 = addDemandReaction(U5, 'cpd00002[c]');
        T5 = addDemandReaction(T5, 'cpd00002[c]');
        U24 = addDemandReaction(U24, 'cpd00002[c]');
        T24 = addDemandReaction(T24, 'cpd00002[c]');
        pa14_dmd = addDemandReaction(pa14, 'cpd00002[c]');
        
        U5 = changeObjective(U5, 'DM_cpd00002[c]',1);
        T5 = changeObjective(T5, 'DM_cpd00002[c]',1);
        U24 = changeObjective(U24, 'DM_cpd00002[c]',1);
        T24 = changeObjective(T24, 'DM_cpd00002[c]',1);
        pa14_dmd = changeObjective(pa14_dmd, 'DM_cpd00002[c]',1);
        
        U5_sol = optimizeCbModel(U5, 'max');
        T5_sol = optimizeCbModel(T5, 'max');
        U24_sol = optimizeCbModel(U24, 'max');
        T24_sol = optimizeCbModel(T24, 'max');
        pa14_sol = optimizeCbModel(pa14_dmd, 'max');
        
        U5.lb(length(U5.rxns)) = 0.1*U5_sol.f;
        T5.lb(length(T5.rxns)) = 0.1*T5_sol.f;
        U24.lb(length(U24.rxns)) = 0.1*U24_sol.f;
        T24.lb(length(T24.rxns)) = 0.1*T24_sol.f;
        pa14_dmd.lb(length(pa14_dmd.rxns)) = 0.1*pa14_sol.f;

        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMin(pa14_dmd, U5, U24, T5, T24);

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
    
        xlswrite('180911_persisterSensitivityAnalysis_MinATP_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
        xlswrite('180911_persisterSensitivityAnalysis_MinATP_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
        xlswrite('180911_persisterSensitivityAnalysis_MinATP_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
        xlswrite('180911_persisterSensitivityAnalysis_MinATP_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

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

  xlswrite('180911_persisterSensitivityAnalysis_minATP_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
  xlswrite('180911_persisterSensitivityAnalysis_minATP_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
  xlswrite('180911_persisterSensitivityAnalysis_minATP_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
  xlswrite('180911_persisterSensitivityAnalysis_minATP_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
  xlswrite('180911_persisterSensitivityAnalysis_minATP_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');

%% Sensitivity Analysis - minimize ATP production w/o minimizing biomass

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
        
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);
        
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
        
        U5 = addDemandReaction(U5, 'cpd00002[c]');
        T5 = addDemandReaction(T5, 'cpd00002[c]');
        U24 = addDemandReaction(U24, 'cpd00002[c]');
        T24 = addDemandReaction(T24, 'cpd00002[c]');
        pa14_dmd = addDemandReaction(pa14, 'cpd00002[c]');
        
        U5 = changeObjective(U5, 'DM_cpd00002[c]',1);
        T5 = changeObjective(T5, 'DM_cpd00002[c]',1);
        U24 = changeObjective(U24, 'DM_cpd00002[c]',1);
        T24 = changeObjective(T24, 'DM_cpd00002[c]',1);
        pa14_dmd = changeObjective(pa14_dmd, 'DM_cpd00002[c]',1);
        
        U5_sol = optimizeCbModel(U5, 'max');
        T5_sol = optimizeCbModel(T5, 'max');
        U24_sol = optimizeCbModel(U24, 'max');
        T24_sol = optimizeCbModel(T24, 'max');
        pa14_sol = optimizeCbModel(pa14_dmd, 'max');
        
        U5.lb(length(U5.rxns)) = 0.1*U5_sol.f;
        T5.lb(length(T5.rxns)) = 0.1*T5_sol.f;
        U24.lb(length(U24.rxns)) = 0.1*U24_sol.f;
        T24.lb(length(T24.rxns)) = 0.1*T24_sol.f;
        pa14_dmd.lb(length(pa14_dmd.rxns)) = 0.1*pa14_sol.f;

        RNA_pval
        RNA_fc
        
        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMin(pa14_dmd, U5, U24, T5, T24);

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
    
%         xlswrite('180921_persisterSensitivityAnalysis_MinATPwoBio_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinATPwoBio_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinATPwoBio_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
%         xlswrite('180921_persisterSensitivityAnalysis_MinATPwoBio_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

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

  xlswrite('180921_persisterSensitivityAnalysis_minATPwoBio_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
  xlswrite('180921_persisterSensitivityAnalysis_minATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
  xlswrite('180921_persisterSensitivityAnalysis_minATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
  xlswrite('180921_persisterSensitivityAnalysis_minATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
  xlswrite('180921_persisterSensitivityAnalysis_minATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');

  %% Sensitivity Analysis - maximize ATP production w/o minimizing biomass

% pval_range = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
% log2FC_range = [0:0.5:2.5];

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
        
        RNA_pval = pval_range(i);
        RNA_fc = log2FC_range(j);
    
        [U5, U24, T5, T24, untreated_metabs_5_int, persister_metabs_5_int, untreated_metabs_24_int, persister_metabs_24_int] = persisterIntegration(pa14, diffEX_file, RNA_pval, RNA_fc);
        
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
        
        U5 = addDemandReaction(U5, 'cpd00002[c]');
        T5 = addDemandReaction(T5, 'cpd00002[c]');
        U24 = addDemandReaction(U24, 'cpd00002[c]');
        T24 = addDemandReaction(T24, 'cpd00002[c]');
        pa14_dmd = addDemandReaction(pa14, 'cpd00002[c]');
        
        U5 = changeObjective(U5, 'DM_cpd00002[c]',1);
        T5 = changeObjective(T5, 'DM_cpd00002[c]',1);
        U24 = changeObjective(U24, 'DM_cpd00002[c]',1);
        T24 = changeObjective(T24, 'DM_cpd00002[c]',1);
        pa14_dmd = changeObjective(pa14_dmd, 'DM_cpd00002[c]',1);
        
        U5_sol = optimizeCbModel(U5, 'max');
        T5_sol = optimizeCbModel(T5, 'max');
        U24_sol = optimizeCbModel(U24, 'max');
        T24_sol = optimizeCbModel(T24, 'max');
        pa14_sol = optimizeCbModel(pa14_dmd, 'max');
        
        U5.ub(length(U5.rxns)) = 0.1*U5_sol.f;
        T5.ub(length(T5.rxns)) = 0.1*T5_sol.f;
        U24.ub(length(U24.rxns)) = 0.1*U24_sol.f;
        T24.ub(length(T24.rxns)) = 0.1*T24_sol.f;
        pa14_dmd.ub(length(pa14_dmd.rxns)) = 0.1*pa14_sol.f;
        
        RNA_pval
        RNA_fc

        [model_ess_pa14, model_ess_U5, model_ess_U24, model_ess_T5, model_ess_T24, common_ess, unique_treated, unique_untreated, unique_U5, unique_U24, unique_T5, unique_T24] = persisterEssentialityMax(pa14_dmd, U5, U24, T5, T24);

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
    
%         xlswrite('180919_persisterSensitivityAnalysis_MaxATPwoBio_essGenes.xlsx',model_ess_U5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'A1');
%         xlswrite('180919_persisterSensitivityAnalysis_MaxATPwoBio_essGenes.xlsx',model_ess_T5,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'B1');
%         xlswrite('180919_persisterSensitivityAnalysis_MaxATPwoBio_essGenes.xlsx',model_ess_U24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'C1');
%         xlswrite('180919_persisterSensitivityAnalysis_MaxATPwoBio_essGenes.xlsx',model_ess_T24,strcat('pval=',num2str(RNA_pval),', log2FC=',num2str(RNA_fc)),'D1');

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

  xlswrite('180919_persisterSensitivityAnalysis_maxATPwoBio_geneStatusTotal.xlsx',pa14.genes,'Sheet1','A1');
  xlswrite('180919_persisterSensitivityAnalysis_maxATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_U5','Sheet1','B1');
  xlswrite('180919_persisterSensitivityAnalysis_maxATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_T5','Sheet1','C1');
  xlswrite('180919_persisterSensitivityAnalysis_maxATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_U24','Sheet1','D1');
  xlswrite('180919_persisterSensitivityAnalysis_maxATPwoBio_geneStatusTotal.xlsx',geneStatusTotal_T24','Sheet1','E1');
