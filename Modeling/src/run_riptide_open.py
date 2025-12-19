import os
import sys
import cobra
from riptide import read_transcription_file, maxfit
import riptide 
from cobra.flux_analysis.variability import find_essential_reactions, find_essential_genes

# Check for correct arguments
if len(sys.argv) != 4:
    print("Usage: python run_riptide_maxfit.py <counts_file.tsv> <model_file.json> <output_dir>")
    sys.exit(1)

# Load inputs
counts_file = sys.argv[1]
model_file = sys.argv[2]
output_dir = sys.argv[3]
print('Counts file : ' + counts_file)
print('Output dir  : ' + output_dir)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load metabolic model
print("Loading metabolic model...")
model = cobra.io.read_sbml_model(model_file)
model.objective='PA14_Biomass'
model.solver='gurobi'

'''
# Define exchange reaction constraints
lb_5 = {'EX_cpd00027_e','EX_cpd00085_e','EX_cpd00117_e','EX_cpd00035_e','EX_cpd00051_e','EX_cpd00041_e','EX_cpd00084_e','EX_cpd00023_e',
'EX_cpd00033_e','EX_cpd00119_e','EX_cpd00322_e','EX_cpd00107_e','EX_cpd00039_e','EX_cpd00060_e','EX_cpd00129_e','EX_cpd00161_e','EX_cpd00069_e',
'EX_cpd00066_e','EX_cpd00550_e','EX_cpd00054_e','EX_cpd00065_e','EX_cpd00156_e','EX_cpd00644_e','EX_cpd08325_e','EX_cpd00246_e',
'EX_cpd00226_e','EX_cpd00654_e','EX_cpd00184_e','EX_cpd00092_e','EX_cpd00249_e','EX_cpd00438_e','EX_cpd00182_e','EX_cpd00082_e',
'EX_cpd00108_e'}

lb_1000 = {'EX_cpd00971_e',
'EX_cpd00099_e',
'EX_cpd00048_e',
'EX_cpd00205_e',
'EX_cpd00009_e',
'EX_cpd00063_e',
'EX_cpd00254_e',
'EX_cpd00034_e',
'EX_cpd01012_e',
'EX_cpd00531_e',
'EX_cpd00011_e',
'EX_cpd00149_e',
'EX_cpd00058_e',
'EX_cpd00021_e',
'EX_cpd10516_e',
'EX_cpd00030_e',
'EX_cpd11574_e',
'EX_cpd00244_e'}

lb_100 = {'EX_cpd00067_e', 'EX_cpd00001_e'}
lb_18_5 = {'EX_cpd00007_e'}
lb_0_01 = {'EX_cpd00635_e'}
'''

# Reset all exchanges
for rxn in model.exchanges:
    rxn.lower_bound = -1000.0
    rxn.upper_bound = 1000
'''
# Apply constraints
print("Applying exchange constraints...")
for rxn_set, bound in [(lb_5, 5), (lb_1000, 1000), (lb_100, 100), (lb_18_5, 18.5), (lb_0_01, 0.01)]:
    for rxn_id in rxn_set:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -bound
            model.reactions.get_by_id(rxn_id).upper_bound = bound
        else:
            print(f"Warning: {rxn_id} not found in model.")
'''
# Read transcription file
print("Processing transcription data...")
transcription_data = read_transcription_file(counts_file, sep=',')

# Run Riptide maxfit
print("Running Riptide maxfit...")
riptide_object_maxfit = maxfit(transcriptome=transcription_data,model=model)
dir= r"/scratch/jmf7ak/persister/outputs/open/"
name=os.path.splitext(os.path.basename(counts_file))[0]
outFile= 'riptide_results_BIOMASS/' + name
current = dir + outFile
riptide.save_output(riptide_obj=riptide_object_maxfit, path=current)

os.chdir(output_dir)

#get essential reactions and gene
model2=riptide_object_maxfit.model
locals()["essential"] = find_essential_reactions(model2)
locals()["essential_id"] = [x.id for x in locals()["essential"]]
with open('essentialrxns.txt', 'a') as outFile:
    for rxn in locals()["essential_id"]:
        outFile.write(rxn + '\n')

locals()["essential_gene"] = find_essential_genes(riptide_object_maxfit.model)
locals()["essential_gene_id"] = [x.id for x in locals()["essential_gene"]]
with open('essentialgenes.txt', 'a') as outFile:
    for gene in locals()["essential_gene_id"]:
        outFile.write(gene + '\n')
