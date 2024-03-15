import os

PROJECT_NAME = "SCA1_iPS"
CASE_GROUP = "SCA1"
CONTROL_GROUP_NAME = "Normal"
CASE_GROUP_NAME = CASE_GROUP

TIME_POINTS = ["iPS", "EB", "Purkinje"]

PROJECT_DIR = "hSCA1_DESeq2"
GENE_EXP_DIR_NAME = "expression_data"
GENE_EXP_DIR = os.path.join(PROJECT_DIR, GENE_EXP_DIR_NAME)
ANNOTATED_GENE_EXP_DIR = os.path.join(PROJECT_DIR, "%s_annotated" % GENE_EXP_DIR_NAME)
CHANGED_MOLECULES_FILE = os.path.join(ANNOTATED_GENE_EXP_DIR, "changed_molecules.tsv")

IMAD_RESULT_DIR = os.path.join(PROJECT_DIR, "iMAD_results")
IMPACT_RESULT_FILE = os.path.join(IMAD_RESULT_DIR, "impact_nodes.tsv")
NORMAL_NODE_FILE = os.path.join(IMAD_RESULT_DIR, "normal_nodes.tsv")

IMAD_NETWORK_HTML = os.path.join(IMAD_RESULT_DIR, "iMAD-based_network_%s.html" % PROJECT_NAME)
IMPACT_EDGES_FILE = os.path.join(IMAD_RESULT_DIR, "impact_edges.tsv")
IMPACT_DOWNSTREAM_FILE = os.path.join(IMAD_RESULT_DIR, "impact_downstream.tsv")
