function model = xls2model_JAB(fileName,biomassRxnEquation)
% xls2model Writes a model from Excel spreadsheet.
%
% model = xls2model(fileName,metFileName)
%
% INPUT
% fileName      xls spreadsheet, with one 'reactions' and one 'metabolites' tab
%
% 'reactions' tab: Required headers:
%               col 1     Abbreviation    HEX1
%               col 2     Name            Hexokinase
%               col 3     Reaction        1 atp[c] + 1 glc-D[c] --> 1 adp[c] + 1 g6p[c] + 1 h[c]
%               col 4     GPR             b0001
%               col 5     Genes           b0001 (optional: column can be empty)
%               col 6     Protein           AlaS (optional: column can be empty)
%               col 7     Subsystem       Glycolysis
%               col 8     Reversible      0
%               col 9     Lower bound     0
%               col 10    Upper bound     1000
%               col 11    Objective       0    (optional: column can be empty)
%               col 12    Confidence Score 0,1,2,3,4
%               col 13    EC. Number      1.1.1.1
%               col 14    Notes           N/A  (optional: column can be empty)
%               col 15    References      PMID: 1111111  (optional: column can be empty)
%
% 'metabolites' tab: Required headers: (needs to be complete list of metabolites, i.e., if a metabolite appears in multiple compartments it has to be represented in multiple rows. Abbreviations needs to overlap with use in Reaction List
%               col 1     Abbreviation
%               col 2     Name
%               col 3     Formula (neutral)
%               col 4     Formula (charged)
%               col 5     Charge
%               col 6     Compartment
%               col 7     KEGG ID
%               col 8     PubChem ID
%               col 9     ChEBI ID
%               col 10    InChI string
%               col 11    Smiles
%
% OPTIONAL INPUT (may be required for input on unix macines)
% biomassRxnEquation        .xls may have a 255 character limit on each cell, 
%                           so pass the biomass reaction separately if it hits this maximum.        
%
% OUTPUT
% model         COBRA Toolbox model
%
% Ines Thiele   01/02/09
% Richard Que   04/27/10    Modified reading of PubChemID and ChEBIID so that if met 
%                           has multiple IDs, all are passed to model. Confidence Scores
%                           PubChemIDs, and ChEBIIDs, are properly passed as cell arrays. 
% Ronan Fleming 08/17/10    Support for unix
% Jennie Bartell 2014

warning('xls2model IS NOT SUPPORTED BY THE openCOBRA CORE TEAM AND WILL BE MOVED FROM THE CORE IN THE NEAR FUTURE');
warning off

if isunix
    %assumes that one has an xls file with two tabs
    [Numbers, Strings, Raw] = xlsread(fileName,'reactions');
    [MetNumbers, MetStrings, MetRaw] = xlsread(fileName,'metabolites');
    %trim empty row from Numbers and MetNumbers
    %Numbers = Numbers(2:end,:);
    %MetNumbers = MetNumbers(2:end,:);
    
    if isempty(MetStrings)
        error('Save .xls file as Windows 95 version using gnumeric not openoffice!');
    end
    
    nRxns=length(Strings(:,1))-1;
    nMets=length(MetStrings(:,1))-1;
    
    %[Numbers, Strings] = xlsread(fileName,'reactions',['A1:O' nRxns],'basic');
    %[MetNumbers, MetStrings] = xlsread(fileName,'metabolites',['A1:K' nMets],'basic');
else
    %assumes that one has an xls file with two tabs
    [Numbers, Strings, Raw] = xlsread(fileName,'reactions');
    assignin('base','Strings',Strings)
    [MetNumbers, MetStrings, MetRaw] = xlsread(fileName,'metabolites');
    assignin('base','MetStrings',MetStrings)
    assignin('base','MetNumbers',MetNumbers)
    % assumed that first row is header row
    nRxns=length(Strings(:,1))-1;
    nMets=length(MetStrings(:,1))-1;
    %add empty line to 
end

rxn_strings={'Abbreviation' 'Name' 'Reaction' 'GPR' 'Genes' 'Protein'...
            'Subsystems' 'Reversible' 'Lower bound' 'Upper bound'...
            'Objective' 'Confidence Score' 'EC. Number' 'Notes' 'References'};
        
met_strings={'Abbreviation' 'Name' 'Formula (neutral)' 'Formula (Charged)' 'Charge'...
            'Compartment' 'KEGG ID' 'PubChem ID' 'ChEBI ID'...
            'InChI string' 'Smiles'};
        
for ctr=1:length(rxn_strings)
    ind=find(ismember(Strings(1,:),rxn_strings(ctr)));
    if ~isempty(ind)
        indR(ctr)=ind;
    else
        indR(ctr)=0;
    end
end

for ctm=1:length(met_strings)
    ind=find(ismember(MetStrings(1,:),met_strings(ctm)));
    if ~isempty(ind)
        indM(ctm)=ind;
    else
        indM(ctm)=0;
    end
end

rxnAbrList = Strings(2:end,indR(1)); 
rxnNameList = Strings(2:end,indR(2));
rxnList = Strings(2:end,indR(3));
grRuleList = Strings(2:end,indR(4));
assignin('base','grRuleList',grRuleList)
subSystemList = Strings(2:end,indR(7));

if isunix
    for n=1:length(rxnList)
        if length(rxnList{n})==255
            if exist('biomassRxnEquation','var')
                rxnList{n}=biomassRxnEquation;
            else
                error('biomassRxnEquation .xls may have a 255 character limit on each cell, so pass the biomass reaction separately if it hits this maximum.')
            end
        end
    end
end

[r,c] = size(Numbers);
%if c >= 1
if indR(6)~=0
    Protein = Strings(2:end,indR(6));
else
    Protein = [];
end

if indR(8)~=0
    revFlagList = cell2mat(Raw(2:end,indR(8)));
else
    revFlagList = [];
end
%if c >= 2
if indR(9)~=0
    lowerBoundList = cell2mat(Raw(2:end,indR(9)));
else
    lowerBoundList = 1000*ones(length(rxnAbrList),1);
end
%if c >= 3
if indR(10)~=0
    upperBoundList = cell2mat(Raw(2:end,indR(10)));
else
    upperBoundList = 1000*ones(length(rxnAbrList),1);
end
%if c >= 4
if indR(11)~=0
    Objective = cell2mat(Raw(2:end,indR(11)));
else
    Objective = zeros(length(rxnAbrList),1);
end

model = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList,lowerBoundList,upperBoundList,subSystemList,grRuleList);
assignin('base','model',model)

if size(Numbers,2)>=5
%if indR(12)~=0
    ConfidenceScore = Numbers(:,5);
    model.confidenceScores = regexprep(cellstr(num2str(ConfidenceScore)),'NaN| ','');
else
    model.confidenceScores = cell(length(model.rxns),1); %empty cell instead of NaN
end
%if size(Strings,2)>=13
if indR(13)~=0
    model.rxnECNumbers = Strings(2:end,indR(13));
end
%if size(Strings,2)>=14
if indR(14)~=0
    model.rxnNotes = Strings(2:end,indR(14));
end
%if size(Strings,2)>=15
if ~isempty(indR(15))
    model.rxnReferences = Strings(2:end,indR(15));
end

%adding in extra columns/notes from spreadsheet outside of standard model
%format, using column header as variable name in COBRA structure
indEV=find(~ismember(Strings(1,:),rxn_strings));
v=genvarname(Strings(1,indEV));

for vind=1:length(v)
    eval([strcat('model.',v{vind}) '= Strings(2:end,indEV(vind));']);
end

indEVm=find(~ismember(MetStrings(1,:),met_strings));
v=genvarname(MetStrings(1,indEVm));
for vind=1:length(v)
    eval([strcat('model.',v{vind}) '= MetStrings(2:end,indEVm(vind));']);
end

%fill in opt info for metabolites
[nMetNum mMetNum] = size(MetNumbers);
if mMetNum<5
    MetNumbers(:,mMetNum+1:5) = nan(nMetNum,5-mMetNum);
end

if ~isempty(Objective) && length(Objective) == length(model.rxns)
    model.c = (Objective);
end
model.proteins = Protein;

% case 1: all metabolites in List have a compartment assignment

if ~cellfun('isempty',(strfind(MetStrings(2,1),'[')))
    for i = 2 : length(MetStrings(:,1))% assumes that first row is header
        % finds metabolites in model structure
        MetLoc =  strmatch(MetStrings{i,1},model.mets,'exact');
        if ~isempty(MetLoc)
            model.metNames{MetLoc} = MetStrings{i,2};
            model.metFormulasNeutral{MetLoc} = MetStrings{i,3};
         %   model.metFormulas{MetLoc} = char(MetStrings{i,4});
            model.metFormulas{MetLoc} = MetStrings{i,4};
            model.metCompartment{MetLoc} = MetStrings{i,6};
            model.metKEGGID{MetLoc} = MetStrings{i,7};
            if size(MetStrings,2) >= 10
                model.metInChIString{MetLoc} = MetStrings{i,10};
            end
            if size(MetStrings,2) >= 11
                model.metSmiles{MetLoc} = MetStrings{i,11};
            end
            if ~isempty(MetNumbers)
                model.metCharge(MetLoc) = MetNumbers(i-1,1);
                if (~isnan(MetNumbers(i-1,4)))
                    model.metPubChemID(MetLoc) = num2Cell(MetNumbers(i-1,4));
                else
                    model.metPubChemID(MetLoc) = MetStrings(i,8);
                end
                if (~isnan(MetNumbers(i-1,5)))
                    model.metChEBIID(MetLoc) = num2Cell(MetNumbers(i-1,5));
                else
                    model.metChEBIID(MetLoc) = MetStrings(i,9);
                end
            end
        else
            warning(['Metabolite ' MetStrings{i,1} ' not in model']);
        end
        MetLoc=[];
    end
else
    % case 2: all metabolites in List have no compartment assignement
    for i = 2 : length(MetStrings(:,1))% assumes that first row is header
        % finds metabolites in model structure
        % this assumes that the compartment is shown with '[ ]'
        MetLoc =  strmatch(strcat(MetStrings{i,1},'['),model.mets);
        if ~isempty(MetLoc)
            for j = 1 : length(MetLoc)
                model.metNames{MetLoc(j)} = MetStrings{i,2};
                model.metFormulasNeutral{MetLoc(j)} = MetStrings{i,3};
                model.metFormulas{MetLoc(j)} = MetStrings{i,4};
                model.metCompartment{MetLoc(j)} = MetStrings{i,6};
                model.metKEGGID{MetLoc(j)} = MetStrings{i,8};
                if size(MetStrings,2) >= 10
                    model.metInChIString{MetLoc(j)} = MetStrings{i,10};
                end
                if size(MetStrings,2) >= 11
                    model.metSmiles{MetLoc(j)} = MetStrings{i,11};
                end
                if ~isempty(MetNumbers)
                    model.metCharge(MetLoc) = MetNumbers(i-1,1);
                    if (~isnan(MetNumbers(i-1,4)))
                        model.metPubChemID(MetLoc) = num2cell(MetNumbers(i-1,4));
                    else
                        model.metPubChemID(MetLoc) = MetStrings(i,8);
                    end
                    if (~isnan(MetNumbers(i-1,5)))
                        model.metChEBIID(MetLoc) = num2Cell(MetNumbers(i-1,5));
                    else
                        model.metChEBIID(MetLoc) = MetStrings(i,9);
                    end
                end
            end
        else
            warning(['Metabolite ' MetStrings{i,1} ' not in model']);
        end
        MetLoc=[];
    end
end

% %% Verify all vectors are column Vectors
model.lb = columnVector(model.lb);
model.ub = columnVector(model.ub);
%model.rev = columnVector(model.rev);
model.c = columnVector(model.c);
model.b = columnVector(model.b);
model.rxns = columnVector(model.rxns);
model.rxnNames = columnVector(model.rxnNames);
model.mets = columnVector(model.mets);
model.metNames = columnVector(model.metNames);
model.metFormulas = columnVector(model.metFormulas);
% commented out since changing this will require updating all model versions used for testing
% model.metCharges = columnVector(model.metCharge); % all others have plural for vector 
if isfield(model,'metCharge')
    model.metCharge = columnVector(model.metCharge);
end
if isfield(model,'metFormulasNeutral')
    model.metFormulasNeutral = columnVector(model.metFormulasNeutral);
end
model.subSystems = columnVector(model.subSystems);
model.rules = columnVector(model.rules);
model.grRules = columnVector(model.grRules);
model.genes = columnVector(model.genes);
model.confidenceScores = columnVector(model.confidenceScores);
model.rxnECNumbers = columnVector(model.rxnECNumbers);
model.rxnNotes = columnVector(model.rxnNotes);
model.rxnReferences = columnVector(model.rxnReferences);
model.proteins = columnVector(model.proteins);
if isfield(model,'metPubChemID')
    model.metPubChemID = columnVector(model.metPubChemID);
end
if isfield(model,'metChEBIID')
    model.metChEBIID = columnVector(model.metChEBIID);
end
if isfield(model,'metCompartment')
    model.metCompartment = columnVector(model.metCompartment);
end
if isfield(model,'metKEGGID')
    model.metKEGGID = columnVector(model.metKEGGID);
end
if isfield(model,'metInChIString')
    model.metInChIString = columnVector(model.metInChIString);
end
if isfield(model,'metSmiles')
    model.metSmiles = columnVector(model.metSmiles);
end

warning on
