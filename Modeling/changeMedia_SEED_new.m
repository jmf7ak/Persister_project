function modelout = changeMedia_SEED(model, mediaCondition, limEx)
%changeMedia -  Changes the lower and upper bounds of the model's exchange
%reactions to alter the in silico media. Current media offerings are LB
%media, glucose minimal media and SCFM.
%
%INPUTS
% model - COBRA model structure
% mediaCondition - A scalar between one and three to indicate the desired 
% in silico media condition. 1 - returns in silico LB. 2 - returns in 
% silico SCFM. 3 - returns in silico glucose minimal media. 3 - returns in silico SCFM2, 
% 5 - returns open exchange meda (e.g., approximate BHI)
%
%OUTPUT
% modelout - COBRA model structure with modified exchange reaction bounds
%
%Matthew Oberhardt, 1-7-2010
%Jennifer Bartell, 3-27-2013
%Anna Blazier, 9-18-2012
%Anna Blazier, 8-8-2018: add open exchanges media condition

%Finds the indices of the exchange reactions in the model
exchangerxns = [];
for rxn = 1:length(model.rxns);
    exchangerxns(rxn) = strncmp('EX_',model.rxns{rxn},3);
end
exchangeindices = find(exchangerxns); 

modelout = model;

modelout.lb(exchangeindices) = zeros(size(exchangeindices));
modelout.ub(exchangeindices) = 1000*ones(size(exchangeindices));
modelout.lb(findRxnIDs(modelout,{'EX_cpd00007_e'})) = -20;  %-1000; %Limit aerobic growth to 20 mmol/gDW/hr of O2 (Added by Phil, EX_O2 removed from openexchanges)
%% LB
if (mediaCondition == 1)

%Nutrients such as ions that are freely exchanged:      
openexchanges = {
'EX_cpd00001_e'    %H2O
'EX_cpd00009_e'    %Phosphate
'EX_cpd00011_e'    %CO2
'EX_cpd00021_e'    %Fe2+
'EX_cpd00034_e'    %Zn2+
'EX_cpd00048_e'    %Sulfate
'EX_cpd00058_e'    %Cu2+
'EX_cpd00205_e'    %K+
'EX_cpd00254_e'    %Mg
'EX_cpd00971_e'    %Na+
'EX_cpd01012_e'    %Cd2+
'EX_cpd00067_e'    %H+
'EX_cpd00528_e' % added SEED Nitrogen (changed by JAB, 01/14/14, not in older LB media formulations)
};

%Limited nutrients in LB
LBexchanges = {
'EX_cpd00023_e'    %L-Glutamate
'EX_cpd00027_e'    %D-Glucose
'EX_cpd00033_e'    %Glycine
'EX_cpd00035_e'    %L-Alanine
'EX_cpd00039_e'    %L-Lysine
'EX_cpd00041_e'    %L-Aspartate
'EX_cpd00051_e'    %L-Arginine
'EX_cpd00054_e'    %L-Serine
'EX_cpd00060_e'    %L-Methionine
'EX_cpd00065_e'    %Tryptophan
'EX_cpd00066_e'    %L-Phenylalanine
'EX_cpd00069_e'    %L-Tyrosine
'EX_cpd00084_e'    %L-Cysteine
'EX_cpd00107_e'    %L-Leucine
'EX_cpd00119_e'    %L-Histidine
'EX_cpd00129_e'    %L-Proline
'EX_cpd00156_e'    %L-Valine
'EX_cpd00161_e'    %L-Threonine
'EX_cpd00305_e'    %Thiamin
'EX_cpd00322_e'    %L-Isoleucine
'EX_cpd00092_e'    %Uracil
'EX_cpd00307_e'    %Cytosine
'EX_cpd03091_e'    %5'-Deoxyadenosine
};

%changes the lower bound of openexchanges to -1000, the lower bound of
%LBexchanges to -10 and the upper bound of LBexchanges to 10.  Also
%changes the upper bound of the glucose exchange reaction to 0.
modelout.lb(find(ismember(model.rxns,openexchanges))) = -1000*ones(size(openexchanges));
modelout.lb(find(ismember(model.rxns,LBexchanges))) = -10*ones(size(LBexchanges));
modelout.ub(find(ismember(model.rxns,LBexchanges))) = 10*ones(size(LBexchanges));
modelout.ub(find(ismember(model.rxns,'EX_cpd00027_e'))) = 0;

%% SCFM
elseif (mediaCondition == 2)

%Nutrients in SCFM
openexchanges = {
'EX_cpd00001_e'    %H2O
'EX_cpd00009_e'    %Phosphate
'EX_cpd00011_e'    %CO2
'EX_cpd00021_e'    %Fe2+
'EX_cpd00023_e'    %L-Glutamate
'EX_cpd00027_e'    %D-Glucose
'EX_cpd00033_e'    %Glycine
'EX_cpd00035_e'    %L-Alanine
'EX_cpd00039_e'    %L-Lysine
'EX_cpd00041_e'    %L-Aspartate
'EX_cpd00048_e'    %Sulfate
'EX_cpd00051_e'    %L-Arginine
'EX_cpd00054_e'    %L-Serine
'EX_cpd00060_e'    %L-Methionine
'EX_cpd00064_e'    %Ornithine
'EX_cpd00065_e'    %Tryptophan
'EX_cpd00066_e'    %L-Phenylalanine
'EX_cpd00067_e'    %H+
'EX_cpd00069_e'    %L-Tyrosine
'EX_cpd00084_e'    %L-Cysteine
'EX_cpd00107_e'    %L-Leucine
'EX_cpd00119_e'    %L-Histidine
'EX_cpd00129_e'    %L-Proline
'EX_cpd00156_e'    %L-Valine
'EX_cpd00221_e'    %D-Lactate
'EX_cpd00161_e'    %L-Threonine
'EX_cpd00205_e'    %K+
'EX_cpd00209_e'    %Nitrate
'EX_cpd00254_e'    %Mg
'EX_cpd00322_e'    %L-Isoleucine
'EX_cpd00971_e'    %Na+
'EX_cpd00013_e'    %NH3

};

ub_check={
'EX_cpd00305_e'    %Thiamin
'EX_cpd00092_e'    %Uracil
'EX_cpd00307_e'    %Cytosine
'EX_cpd03091_e'};    %5'-Deoxyadenosine

%changes the lower bound for all available nutrients to -10.
modelout.lb(find(ismember(model.rxns,openexchanges))) = -10*ones(size(openexchanges));
%modelout.ub(find(ismember(model.rxns,[openexchanges;ub_check]))) = 10*ones(size([openexchanges;ub_check]));
%modelout.ub(find(ismember(model.rxns,[openexchanges]))) = 10*ones(size([openexchanges]));
%modelout.ub(find(ismember(model.rxns,'EX_cpd00027_e'))) = 0;
%% Minimal Media
elseif (mediaCondition == 3)

%Nutrients such as ions that are freely exchanged:    
openexchanges = {     
 'EX_cpd00001_e' %H2O
 'EX_cpd00009_e' %Phosphate
 'EX_cpd00011_e' %CO2
 'EX_cpd00021_e' %Fe2+
 'EX_cpd00030_e' %Mn2+
 'EX_cpd00034_e' %Zn2+
 'EX_cpd00048_e' %Sulfate
 'EX_cpd00058_e' %Cu2+
 'EX_cpd00067_e' %H+
 'EX_cpd00149_e' %Co2+
 'EX_cpd00205_e' %K+
 'EX_cpd00254_e' %Mg
 'EX_cpd00528_e' %SEED Nitrogen (changed by JAB, 01/14/14) %'EX_cJB00102_e' %Nitrogen
 'EX_cpd00971_e' %Na+
 'EX_cpd00013_e' %NH3
 'EX_cpd01012_e' %Cd2+
 'EX_cpd10516_e' %fe3
 'EX_cpd00244_e' %Ni2+ ADDED BY PHIL
 %'EX_cpd00209_e' %Nitrate Added by JAB 04/05/14 - checking anaerobic pathways
};

% openexchanges = {''}

%Minimal media carbon source exchange reaction:
 %limitedexchanges = {  
 %'EX_cpd00027_e'}; %glc
%limitedexchanges = {''};

limitedexchanges = limEx;

%changes the lower bound of openexchanges to -1000 and the lower bound of
%the limitedexchanges to -10.  Also changes the upper bound of
%limitedexchanges to 0.
modelout.lb(find(ismember(model.rxns,openexchanges))) = -1000*ones(size(openexchanges));
modelout.lb(find(ismember(model.rxns,limitedexchanges))) = -10;%*ones(size(limitedexchanges));
%modelout.ub(find(ismember(model.rxns,limitedexchanges))) = zeros(size(limitedexchanges));

%% SCFM2
elseif (mediaCondition == 4)

%Nutrients in SCFM
openexchanges = {
'EX_cpd00001_e'    %H2O
'EX_cpd00009_e'    %Phosphate
'EX_cpd00011_e'    %CO2
'EX_cpd00021_e'    %Fe2+
'EX_cpd00023_e'    %L-Glutamate
'EX_cpd00027_e'    %D-Glucose
'EX_cpd00033_e'    %Glycine
'EX_cpd00035_e'    %L-Alanine
'EX_cpd00039_e'    %L-Lysine
'EX_cpd00048_e'    %Sulfate
'EX_cpd00051_e'    %L-Arginine
'EX_cpd00054_e'    %L-Serine
'EX_cpd00060_e'    %L-Methionine
'EX_cpd00064_e'    %Ornithine
'EX_cpd00065_e'    %Tryptophan
'EX_cpd00066_e'    %L-Phenylalanine
'EX_cpd00067_e'    %H+
'EX_cpd00084_e'    %L-Cysteine
'EX_cpd00107_e'    %L-Leucine
'EX_cpd00119_e'    %L-Histidine
'EX_cpd00129_e'    %L-Proline
'EX_cpd00156_e'    %L-Valine
'EX_cpd00221_e'    %D-Lactate
'EX_cpd00161_e'    %L-Threonine
'EX_cpd00205_e'    %K+
'EX_cpd00209_e'    %Nitrate
'EX_cpd00254_e'    %Mg
'EX_cpd00322_e'    %L-Isoleucine
'EX_cpd00971_e'    %Na+
'EX_cpd00013_e'    %NH3

};

ub_check={
'EX_cpd00305_e'    %Thiamin
'EX_cpd00092_e'    %Uracil
'EX_cpd00307_e'    %Cytosine
'EX_cpd03091_e'};    %5'-Deoxyadenosine

%changes the lower bound for all available nutrients to -10.
modelout.lb(find(ismember(model.rxns,openexchanges))) = -10*ones(size(openexchanges));
%modelout.ub(find(ismember(model.rxns,[openexchanges;ub_check]))) = 10*ones(size([openexchanges;ub_check]));
%modelout.ub(find(ismember(model.rxns,[openexchanges]))) = 10*ones(size([openexchanges]));
%modelout.ub(find(ismember(model.rxns,'EX_cpd00027_e'))) = 0;

%% Open exchange (e.g., approximate BHI)
elseif (mediaCondition == 5)
    
modelout.lb(exchangeindices) = -1000*ones(size(exchangeindices));

end

modelout;



end
