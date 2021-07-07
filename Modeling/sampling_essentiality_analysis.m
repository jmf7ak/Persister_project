clear all
close all
clc

pa14 = readCbModel('iPau21');
pa14 = creategrRulesField(pa14);
%% Basal biomass genes of interest

% Pull out the genes associated with the reactions that have a greater
% likelihood of carrying flux in treated conditions

[ndata,sampling_Bio_median_rxns_T5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_median_posthoc_sig_T5_rxns.xlsx');

sampling_Bio_median_genes_T5 = [];

for i = 1:length(sampling_Bio_median_rxns_T5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_median_rxns_T5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_median_genes_idx = find(rxn_genes == 1);
    sampling_Bio_median_genes_T5 = [sampling_Bio_median_genes_T5;pa14.genes(sampling_Bio_median_genes_idx)];
end

sampling_Bio_median_genes_T5 = unique(sampling_Bio_median_genes_T5);

sampling_Bio_median_genes_T5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_median_genes_T5,pa14.genes(i))) ~= 0
        sampling_Bio_median_genes_T5_write(i) = {'yes'};
    end
end
    
[ndata,sampling_Bio_median_rxns_T24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_median_posthoc_sig_T24_rxns.xlsx');

sampling_Bio_median_genes_T24 = [];

for i = 1:length(sampling_Bio_median_rxns_T24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_median_rxns_T24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_median_genes_idx = find(rxn_genes == 1);
    sampling_Bio_median_genes_T24 = [sampling_Bio_median_genes_T24;pa14.genes(sampling_Bio_median_genes_idx)];
end

sampling_Bio_median_genes_T24 = unique(sampling_Bio_median_genes_T24);

sampling_Bio_median_genes_T24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_median_genes_T24,pa14.genes(i))) ~= 0
        sampling_Bio_median_genes_T24_write(i) = {'yes'};
    end
end

[ndata,sampling_Bio_median_rxns_T,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_median_posthoc_sig_T_rxns.xlsx');

sampling_Bio_median_genes_T = [];

for i = 1:length(sampling_Bio_median_rxns_T)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_median_rxns_T(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_median_genes_idx = find(rxn_genes == 1);
    sampling_Bio_median_genes_T = [sampling_Bio_median_genes_T;pa14.genes(sampling_Bio_median_genes_idx)];
end

sampling_Bio_median_genes_T = unique(sampling_Bio_median_genes_T);

sampling_Bio_median_genes_T_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_median_genes_T,pa14.genes(i))) ~= 0
        sampling_Bio_median_genes_T_write(i) = {'yes'};
    end
end

% Pull out the genes associated with the reactions that are highly
% correlated with essential reactions

[ndata,sampling_Bio_correlated_rxns_T5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_corr_rxns_T5.xlsx');

sampling_Bio_correlated_genes_T5 = [];

for i = 1:length(sampling_Bio_correlated_rxns_T5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_correlated_rxns_T5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_correlated_genes_idx = find(rxn_genes == 1);
    sampling_Bio_correlated_genes_T5 = [sampling_Bio_correlated_genes_T5;pa14.genes(sampling_Bio_correlated_genes_idx)];
end

sampling_Bio_correlated_genes_T5 = unique(sampling_Bio_correlated_genes_T5);

sampling_Bio_correlated_genes_T5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_correlated_genes_T5,pa14.genes(i))) ~= 0
        sampling_Bio_correlated_genes_T5_write(i) = {'yes'};
    end
end

[ndata,sampling_Bio_correlated_rxns_T24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_corr_rxns_T24.xlsx');

sampling_Bio_correlated_genes_T24 = [];

for i = 1:length(sampling_Bio_correlated_rxns_T24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_correlated_rxns_T24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_correlated_genes_idx = find(rxn_genes == 1);
    sampling_Bio_correlated_genes_T24 = [sampling_Bio_correlated_genes_T24;pa14.genes(sampling_Bio_correlated_genes_idx)];
end

sampling_Bio_correlated_genes_T24 = unique(sampling_Bio_correlated_genes_T24);

sampling_Bio_correlated_genes_T24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_correlated_genes_T24,pa14.genes(i))) ~= 0
        sampling_Bio_correlated_genes_T24_write(i) = {'yes'};
    end
end

[ndata,sampling_Bio_correlated_rxns_U5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_corr_rxns_U5.xlsx');

sampling_Bio_correlated_genes_U5 = [];

for i = 1:length(sampling_Bio_correlated_rxns_U5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_correlated_rxns_U5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_correlated_genes_idx = find(rxn_genes == 1);
    sampling_Bio_correlated_genes_U5 = [sampling_Bio_correlated_genes_U5;pa14.genes(sampling_Bio_correlated_genes_idx)];
end

sampling_Bio_correlated_genes_U5 = unique(sampling_Bio_correlated_genes_U5);

sampling_Bio_correlated_genes_U5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_correlated_genes_U5,pa14.genes(i))) ~= 0
        sampling_Bio_correlated_genes_U5_write(i) = {'yes'};
    end
end

[ndata,sampling_Bio_correlated_rxns_U24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/Bio_corr_rxns_U24.xlsx');

sampling_Bio_correlated_genes_U24 = [];

for i = 1:length(sampling_Bio_correlated_rxns_U24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_Bio_correlated_rxns_U24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_Bio_correlated_genes_idx = find(rxn_genes == 1);
    sampling_Bio_correlated_genes_U24 = [sampling_Bio_correlated_genes_U24;pa14.genes(sampling_Bio_correlated_genes_idx)];
end

sampling_Bio_correlated_genes_U24 = unique(sampling_Bio_correlated_genes_U24);

sampling_Bio_correlated_genes_U24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_Bio_correlated_genes_U24,pa14.genes(i))) ~= 0
        sampling_Bio_correlated_genes_U24_write(i) = {'yes'};
    end
end













% Essentiality anlysis
[ndata,essentiality_Bio_T5_genes,rdata] = xlsread('/Users/papinlab/Dropbox/Lab/Projects/Persisters/Model/essentiality_Bio_T5_genes.xlsx');

essentiality_Bio_T5_rxns = [];
essentiality_Bio_T5_rxnNames = [];
essentiality_Bio_T5_subsys = [];

for i = 1:length(essentiality_Bio_T5_genes)
    gene_idx = find(strcmp(pa14.genes,essentiality_Bio_T5_genes(i)));
    gene_rxns = pa14.rxnGeneMat(:,gene_idx);
    essentiality_Bio_T5_rxns_idx = find(gene_rxns == 1);
    essentiality_Bio_T5_rxns_idx = essentiality_Bio_T5_rxns_idx(1);
    essentiality_Bio_T5_rxns = [essentiality_Bio_T5_rxns;pa14.rxns(essentiality_Bio_T5_rxns_idx)];
    essentiality_Bio_T5_rxnNames = [essentiality_Bio_T5_rxnNames;pa14.rxnNames(essentiality_Bio_T5_rxns_idx)];
    essentiality_Bio_T5_subsys = [essentiality_Bio_T5_subsys;pa14.KEGGSubsys(essentiality_Bio_T5_rxns_idx)];
end

[~,essentiality_Bio_T24_genes,rdata] = xlsread('/Users/papinlab/Dropbox/Lab/Projects/Persisters/Model/essentiality_Bio_T24_genes.xlsx');

essentiality_Bio_T24_rxns = [];
essentiality_Bio_T24_rxnNames = [];
essentiality_Bio_T24_subsys = [];

for i = 1:length(essentiality_Bio_T24_genes)
    gene_idx = find(strcmp(pa14.genes,essentiality_Bio_T24_genes(i)));
    gene_rxns = pa14.rxnGeneMat(:,gene_idx);
    essentiality_Bio_T24_rxns_idx = find(gene_rxns == 1);
    essentiality_Bio_T24_rxns_idx = essentiality_Bio_T24_rxns_idx(1);
    essentiality_Bio_T24_rxns = [essentiality_Bio_T24_rxns;pa14.rxns(essentiality_Bio_T24_rxns_idx)];
    essentiality_Bio_T24_rxnNames = [essentiality_Bio_T24_rxnNames;pa14.rxnNames(essentiality_Bio_T24_rxns_idx)];
    essentiality_Bio_T24_subsys = [essentiality_Bio_T24_subsys;pa14.KEGGSubsys(essentiality_Bio_T24_rxns_idx)];
end

%% Basal ATP genes of interest

% Pull out the genes associated with the reactions that have a greater
% likelihood of carrying flux in treated conditions

[ndata,sampling_ATP_median_rxns_T5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_median_posthoc_sig_T5_rxns.xlsx');

sampling_ATP_median_genes_T5 = [];

for i = 1:length(sampling_ATP_median_rxns_T5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_median_rxns_T5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_median_genes_idx = find(rxn_genes == 1);
    sampling_ATP_median_genes_T5 = [sampling_ATP_median_genes_T5;pa14.genes(sampling_ATP_median_genes_idx)];
end

sampling_ATP_median_genes_T5 = unique(sampling_ATP_median_genes_T5);

sampling_ATP_median_genes_T5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_median_genes_T5,pa14.genes(i))) ~= 0
        sampling_ATP_median_genes_T5_write(i) = {'yes'};
    end
end
    
[ndata,sampling_ATP_median_rxns_T24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_median_posthoc_sig_T24_rxns.xlsx');

sampling_ATP_median_genes_T24 = [];

for i = 1:length(sampling_ATP_median_rxns_T24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_median_rxns_T24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_median_genes_idx = find(rxn_genes == 1);
    sampling_ATP_median_genes_T24 = [sampling_ATP_median_genes_T24;pa14.genes(sampling_ATP_median_genes_idx)];
end

sampling_ATP_median_genes_T24 = unique(sampling_ATP_median_genes_T24);

sampling_ATP_median_genes_T24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_median_genes_T24,pa14.genes(i))) ~= 0
        sampling_ATP_median_genes_T24_write(i) = {'yes'};
    end
end

[ndata,sampling_ATP_median_rxns_T,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_median_posthoc_sig_T_rxns.xlsx');

sampling_ATP_median_genes_T = [];

for i = 1:length(sampling_ATP_median_rxns_T)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_median_rxns_T(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_median_genes_idx = find(rxn_genes == 1);
    sampling_ATP_median_genes_T = [sampling_ATP_median_genes_T;pa14.genes(sampling_ATP_median_genes_idx)];
end

sampling_ATP_median_genes_T = unique(sampling_ATP_median_genes_T);

sampling_ATP_median_genes_T_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_median_genes_T,pa14.genes(i))) ~= 0
        sampling_ATP_median_genes_T_write(i) = {'yes'};
    end
end

% Pull out the genes associated with the reactions that are highly
% correlated with essential reactions

[ndata,sampling_ATP_correlated_rxns_T5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_corr_rxns_T5.xlsx');

sampling_ATP_correlated_genes_T5 = [];

for i = 1:length(sampling_ATP_correlated_rxns_T5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_correlated_rxns_T5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_correlated_genes_idx = find(rxn_genes == 1);
    sampling_ATP_correlated_genes_T5 = [sampling_ATP_correlated_genes_T5;pa14.genes(sampling_ATP_correlated_genes_idx)];
end

sampling_ATP_correlated_genes_T5 = unique(sampling_ATP_correlated_genes_T5);

sampling_ATP_correlated_genes_T5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_correlated_genes_T5,pa14.genes(i))) ~= 0
        sampling_ATP_correlated_genes_T5_write(i) = {'yes'};
    end
end

[ndata,sampling_ATP_correlated_rxns_T24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_corr_rxns_T24.xlsx');

sampling_ATP_correlated_genes_T24 = [];

for i = 1:length(sampling_ATP_correlated_rxns_T24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_correlated_rxns_T24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_correlated_genes_idx = find(rxn_genes == 1);
    sampling_ATP_correlated_genes_T24 = [sampling_ATP_correlated_genes_T24;pa14.genes(sampling_ATP_correlated_genes_idx)];
end

sampling_ATP_correlated_genes_T24 = unique(sampling_ATP_correlated_genes_T24);

sampling_ATP_correlated_genes_T24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_correlated_genes_T24,pa14.genes(i))) ~= 0
        sampling_ATP_correlated_genes_T24_write(i) = {'yes'};
    end
end

[ndata,sampling_ATP_correlated_rxns_U5,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_corr_rxns_U5.xlsx');

sampling_ATP_correlated_genes_U5 = [];

for i = 1:length(sampling_ATP_correlated_rxns_U5)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_correlated_rxns_U5(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_correlated_genes_idx = find(rxn_genes == 1);
    sampling_ATP_correlated_genes_U5 = [sampling_ATP_correlated_genes_U5;pa14.genes(sampling_ATP_correlated_genes_idx)];
end

sampling_ATP_correlated_genes_U5 = unique(sampling_ATP_correlated_genes_U5);

sampling_ATP_correlated_genes_U5_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_correlated_genes_U5,pa14.genes(i))) ~= 0
        sampling_ATP_correlated_genes_U5_write(i) = {'yes'};
    end
end

[ndata,sampling_ATP_correlated_rxns_U24,rdata] = xlsread('/Users/papinlab/Documents/Lab/Persister/Model/ATP_corr_rxns_U24.xlsx');

sampling_ATP_correlated_genes_U24 = [];

for i = 1:length(sampling_ATP_correlated_rxns_U24)
    rxn_idx = find(strcmp(pa14.rxns,sampling_ATP_correlated_rxns_U24(i)));
    rxn_genes = pa14.rxnGeneMat(rxn_idx,:);
    sampling_ATP_correlated_genes_idx = find(rxn_genes == 1);
    sampling_ATP_correlated_genes_U24 = [sampling_ATP_correlated_genes_U24;pa14.genes(sampling_ATP_correlated_genes_idx)];
end

sampling_ATP_correlated_genes_U24 = unique(sampling_ATP_correlated_genes_U24);

sampling_ATP_correlated_genes_U24_write = cell(length(pa14.genes),1);

for i = 1:length(pa14.genes)
    if find(strcmp(sampling_ATP_correlated_genes_U24,pa14.genes(i))) ~= 0
        sampling_ATP_correlated_genes_U24_write(i) = {'yes'};
    end
end

%% Overlap genes of interest

median_genes = intersect(sampling_Bio_median_genes_T5,sampling_ATP_median_genes);

sampling_median_rxns = [];

for i = 1:length(median_genes)
    gene_idx = find(strcmp(pa14.genes,median_genes(i)));
    gene_rxns = pa14.rxnGeneMat(:,gene_idx);
    sampling_median_rxn_idx = find(gene_rxns == 1);
    sampling_median_rxns = [sampling_median_rxns;pa14.rxns(sampling_median_rxn_idx)];
end
 



rxns = [];
rxnNames = [];
subsys = [];

for i = 1:length(pa14.genes)
    gene_idx = i;
    gene_rxns = pa14.rxnGeneMat(:,gene_idx);
    rxns_idx = find(gene_rxns == 1);
    rxns_idx = rxns_idx(1);
    rxns = [rxns;pa14.rxns(rxns_idx)];
    rxnNames = [rxnNames;pa14.rxnNames(rxns_idx)];
    subsys = [subsys;pa14.KEGGSubsys(rxns_idx)];
end
