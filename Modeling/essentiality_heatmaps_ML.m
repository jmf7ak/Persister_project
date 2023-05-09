%% heatmap for essentiality
clc; clear;
cd('C:\Users\jmfic\OneDrive\Documents\joe_outputs\BIOMASS')
%genes
Dir = 'C:\Users\jmfic\OneDrive\Documents\joe_outputs\BIOMASS\genes\';
myFiles = dir(fullfile(Dir));
B=cell(47,1);
for k = 1:length(myFiles)
     baseFileName = myFiles(k).name;
     fullFileName = fullfile(Dir, baseFileName);
     B{k}=baseFileName;
end
B(2)=[];
B(1)=[];
%%
%% Getting file names
R=struct;
for i=1:45
    name=B{i};
    fid = fopen(name);
    data = textscan(fid,'%s%s%s');
    fclose(fid);
    R(i).rxns=data{1,1};
    R(i).name=B{i};
end

%%
allgenes=cell(0,0);
for i=1:45
    allgenes=union(allgenes,R(i).rxns);
end 
%%
genemat=zeros(length(allgenes),45);
for j=1:45
    for i=1:numel(allgenes)
        if contains(allgenes(i),R(j).rxns)
            %spot=find(strcmp(R(j).rxns,allgenes(i)));
            genemat(i,j)=0;
        else 
            genemat(i,j)=-1;
        end 
    end
end

%%
for i=1:45
    newname=strsplit(R(i).name,'.');
    S(i).name=newname{1};
end
%%
for i=250:-1:1
    if sum(genemat(i,:))==0
        genemat(i,:)=[];
        allgenes(i)=[];
    end
end

%%
reactions=clustergram(genemat,'ColumnLabels',{S(:).name},'Colormap',flipud(gray));
addTitle(reactions,'Replicate essential genes')


%% 
% 
% 
% heatmap for RXN essentiality
%rxns
clc; clear;
Dir = 'C:\Users\jmfic\OneDrive\Documents\joe_outputs\BIOMASS\rxns\';
myFiles = dir(fullfile(Dir));
B=cell(47,1);
for k = 1:length(myFiles)
     baseFileName = myFiles(k).name;
     fullFileName = fullfile(Dir, baseFileName);
     B{k}=baseFileName;
end
B(2)=[];
B(1)=[];
%%
%% Getting file names
R=struct;
for i=1:45
    name=B{i};
    fid = fopen(name);
    data = textscan(fid,'%s%s%s');
    fclose(fid);
    R(i).rxns=data{1,1};
    R(i).name=B{i};
end

%%
allrxns=cell(0,0);
for i=1:45
    allrxns=union(allrxns,R(i).rxns);
end 
%%
rxnmat=zeros(length(allrxns),45);
for j=1:45
    for i=1:numel(allrxns)
        if contains(allrxns(i),R(j).rxns)
            %spot=find(strcmp(R(j).rxns,allgenes(i)));
            rxnmat(i,j)=0;
        else 
            rxnmat(i,j)=-1;
        end 
    end
end

%%
rxn_name={};
for i=1:45
    newname=strsplit(R(i).name,'_');
    S(i).name=newname{1};
    rxn_name{i}=newname{1};
end
%%
for i=458:-1:1
    if sum(rxnmat(i,:))==0
        rxnmat(i,:)=[];
        allrxns(i)=[];
    end
end
%%
rxn_name(29:-1:24)=[];
rxnmat(:,29:-1:24)=[];
%%
for i=1:8
    rxn_name(i)={'U0'};
end
%%
reactions=clustergram(rxnmat,'ColumnLabels',rxn_name,'Colormap',flipud(gray));
addTitle(reactions,'Replicate essential reactions')
%%
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 14)





%% get reaction names from model
initCobraToolbox()
%%
pa14 = readCbModel('iPau21.xml');
%%
sol = optimizeCbModel(pa14);
%%
rxns={'rxn00178','rxn01730','rxn00872','rxn10074','rxn13837','EX_cpd00379_e'};

for i=1:187
    spot=find(strcmp(pa14.rxns,allrxns(i)));
    rxn_name(i)=pa14.rxnNames(spot);
end

%%
%%
rxns={'rxn13822','rxn13705','rxn13682','rxn00216','rxn05147','EX_cpd03091_e','EX_cpd00027_e'};

for i=1:7
    spot=find(strcmp(pa14.rxns,rxns(i)));
    rxn_name(i)=pa14.rxnNames(spot);
end




%% MACHINE LEARNING
clc; clear;
%%
R=struct;
Dir = 'C:\Users\jmfic\OneDrive\Documents\joe_outputs\BIOMASS\samples';
j=1;
myFiles = dir(fullfile(Dir));
B=cell(45,1);
for k = 1:length(myFiles)
     baseFileName = myFiles(k).name;
     fullFileName = fullfile(Dir, baseFileName);
     B{k}=baseFileName;
end
B(2)=[];
B(1)=[];
%% R struct has all the flux sampling in it
R=struct;
Dir = 'C:\Users\jmfic\OneDrive\Documents\joe_outputs\BIOMASS\samples\';
for i=1:45
    myDir=strcat(Dir,B{i});
    cd(myDir)
    [fluxes,rxns]=tsvread('flux_samples.tsv');
    R(i).rxns=rxns(1:end-1);
    R(i).flux=fluxes(2:end,2:end);
    R(i).name=B{i};
end

%% get a matrix of all reactions-- 486
all_reactions=R(1).rxns(1:end);
for i=1:45
    all_reactions=union(all_reactions,R(i).rxns(1:end));
end
%%
clc;
MLinput=zeros(900,469);
for i=1:469
    for j=1:45
        location=find(strcmp(R(j).rxns(1:end),char(all_reactions(i))));
        if ~isempty(location)
        MLinput(20*j-19:20*j,i)=R(j).flux(1:20,location);
        else
        MLinput(20*j-19:20*j,i)=0;
        end
    end
end

%%
ML_output=cell(900,1);
ML_output(1:160)={'U0'};
ML_output(161:300)={'P24'};
ML_output(301:460)={'P5'};
ML_output(461:580)={'U0'};
ML_output(581:740)={'U24'};
ML_output(741:900)={'U5'};

%%
impOOB = oobPermutedPredictorImportance(biomass_forest_2.ClassificationEnsemble);

%%
figure
bar(impOOB)
title('Unbiased Predictor Importance Estimates')
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel =trainedModel_persister_biomass.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

%%
pa14 = readCbModel('iPau21.xml');
all_reactions_names=cell(469,1);
for i=1:469
    spot=find(strcmp(pa14.rxns,all_reactions(i)));
    all_reactions_names(i)=pa14.rxnNames(spot);
end

%%
view(biomass_tree_2.ClassificationTree,'mode','graph')
 

%% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(MLinput);
g=gscatter(score(:,1),score(:,2),ML_output,'rkgbmy','.',12);
for i = 1:length(g)
    g(i).MarkerSize = 10;
end

%% bray curtis
dis=f_braycurtis(MLinput');

%% nmds
clc;
Y = f_nmds(dis,2);
%%
gscatter(Y.scores(:,1),Y.scores(:,2),ML_output,'rkgbmy','.',12)
xlabel('NMDS1')
ylabel('NMDS2')

%% delete U0
for i=580:-1:461
    MLinput(i,:)=[];
    ML_output(i,:)=[];
end 

%% TSNE
Y=tsne(MLinput,'Algorithm','exact','Distance','euclidean');
%%
g=gscatter(Y(:,1),Y(:,2),ML_output);

for i = 1:length(g)
    g(i).MarkerSize = 10;
end

xlabel('tSNE1')
ylabel('tSNE2')