import os

VERSION_DATE = "20240126"

TIME_POINTS = ["iPS", "EB", "Purkinje"]

P_THRESHOLD = 0.05
LFC_THRESHOLD = 0.5
IMPACT_RATIO_THRESHOLD = {'iPS': 0, 'EB': 0, 'Purkinje': 0}
IMPACT_P_THRESHOLD = 0.05

SAMPLE_INFO_FILE = "SCA1_meta-analysis_dataset.xlsx"
SAMPLE_INFO_SHEETNAME = "samples_used"

COUNT_DIR_NAME = "featureCounts"
COUNT_DATA_FILE = os.path.join(COUNT_DIR_NAME, "counts.txt")
COUNT_SUMMARY_FILE = os.path.join(COUNT_DIR_NAME, "counts.txt.summary")

ATTRIBUTE_DIR_NAME = os.path.join("network_attributes")
ANNOTATION_FILE_NAME = "uniprot_annotation.tsv" #"gene_id_map_mod.txt"
ANNOTATION_FILE = os.path.join(ATTRIBUTE_DIR_NAME, ANNOTATION_FILE_NAME)
NODE_ATTRIBUTE_FILE_NAME = "SCA1_meta_network_attributes.xlsx"
NODE_ATTRIBUTE_FILE = os.path.join(ATTRIBUTE_DIR_NAME, NODE_ATTRIBUTE_FILE_NAME) # time_points, node_trace, pos_xy
UNIPROT_ANNOTATION_FILE = os.path.join(ATTRIBUTE_DIR_NAME, "uniprot_annotation.tsv")

MOLECULES_DIR_NAME = "molecules"
GENE_ID_MAP_FILE_NAME = "gene_annotation_m2h_all.tsv"
EDGELIST_FILE_NAME = "edge_sca1_integrated.tsv"
MOLECULES_DIR = os.path.join(".", MOLECULES_DIR_NAME)
GENE_ID_MAP_FILE = os.path.join(MOLECULES_DIR, GENE_ID_MAP_FILE_NAME)
ALL_MOLECULES_FILE = os.path.join(MOLECULES_DIR, "all_molecules.tsv")
EDGELIST_FILE = os.path.join(MOLECULES_DIR, EDGELIST_FILE_NAME)

# parameters for gene expression data file
DATA_TYPE = "mRNA_expression" ## mRNA_expression / phospho_proteome
GENE_ID_COLNAME = "Geneid"
UNIPROT_ACC_COLNAME = "uniprot_accession"
GENE_SYMBOL_COLNAME = "symbol"
EXPRESSION_VALUE_COLNAME = "ratio"
PVALUES_COLNAME = "pValue"
## optional:
QVALUES_COLNAME = "qValue"
PROTEIN_NAMES_COLNAME = "symbol"
CASE_VALUE_COLNAME = "mean_case"
CONTROL_VALUE_COLNAME = "mean_control"
