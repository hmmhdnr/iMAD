import os
import sys
import math
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm
from importlib import import_module

from scipy import stats
import statsmodels.stats.multitest as multi

from DyMolNet import load_input_data as dyInput

import common_params as common
import params_hSCA1_DESeq2_LFC05 as param

# init
os.makedirs(param.IMAD_RESULT_DIR, exist_ok=True)

args = {
    "project_name": param.PROJECT_NAME,
    "group_name": param.CASE_GROUP_NAME,
    "data_type": common.DATA_TYPE,
    "time_points": param.TIME_POINTS,
    "data_dir": param.ANNOTATED_GENE_EXP_DIR,
    "id_col": common.GENE_ID_COLNAME,
    "uniprot_col": common.UNIPROT_ACC_COLNAME,
    "symbol_col": common.GENE_SYMBOL_COLNAME,
    "values_col": common.EXPRESSION_VALUE_COLNAME,
    "pvalues_col": common.PVALUES_COLNAME
}
df = dyInput.load_data(**args)

if "uniprot_accession" not in df.columns:
    print("add uniprot accession to input data")
    gene_id_path = common.ANNOTATION_FILE
    id_table = pd.read_table(gene_id_path, sep = "\t").set_index("Gene_ID")
    df["uniprot_accession"] = id_table["uniprot_accession"]
    df["gene_names"] = id_table["symbol"]

if "gene_names" not in df.columns:
    df["gene_names"] = df["gene_symbol"]

os.makedirs(common.MOLECULES_DIR_NAME, exist_ok=True)
molecules, changed, changed_molecules = dyInput.get_molecule_list_lfc(
    df, pvalues_col = "p_value", lfc_col = "value", lfc_threshold = common.LFC_THRESHOLD,
    all_molecule_path = common.ALL_MOLECULES_FILE,
    changed_molecule_path = param.CHANGED_MOLECULES_FILE
)
print("All the molecule associated with UniProt accession number = %d" % len(molecules))
print("number of changed nodes through analysis = %d" % changed_molecules.shape[0])

print("<total candidates>")
for time_point in param.TIME_POINTS:
    _d = changed.loc[changed.time_point == time_point, :]
    print("number of changed molecules at %s = %d" % (time_point, _d.shape[0]))

print("<unique molecules>")
for time_point in param.TIME_POINTS:
    _d = changed_molecules.loc[changed_molecules.time_point == time_point, :]
    print("number of changed molecules at %s = %d" % (time_point, _d.shape[0]))

edgelist = pd.read_table(
    common.EDGELIST_FILE, sep = "\t", index_col = False, header = 0,
    names = ["source", "target"]
)

for i in range(edgelist.shape[0]):
    if edgelist.source[i] > edgelist.target[i]:
        tmp = edgelist.source[i]
        edgelist.source[i] = edgelist.target[i]
        edgelist.target[i] = tmp

edge_id = edgelist.source.str.cat(edgelist.target, sep = "_")
edgelist.loc[:, "edge_id"] = edge_id
edgelist.set_index("edge_id", inplace = True)

def downstream_nodes_direct(G, node):
    p = nx.single_source_shortest_path(G, node, 1)
    n = list(set(list(p.keys())).difference([node]))
    return n

G = nx.from_pandas_edgelist(edgelist)
G = G.to_undirected()
molecules_in_G = [x for x in molecules if x in G]

impact = pd.DataFrame(
    index = [],
    columns = ['uniprot_accession', 'time_point', 'Ratio', 'changed', 'total', 'changed_nodes', 'other_nodes']
)
other_nodes = pd.DataFrame(
    index = [], columns = ['node_id', 'uniprot_accession', 'time_point']
)

for i in tqdm(range(len(param.TIME_POINTS) - 1)):
    tp = param.TIME_POINTS[i]
    tp_next = param.TIME_POINTS[i+1]
    downstream = []
    
    tmp = changed_molecules.loc[changed_molecules.time_point == tp, :]
    print("%d changed proteins at %s" % (tmp.shape[0], tp))
    _candidate_molecules = sorted(list(set(tmp.uniprot_accession)))
    print("%d unique candidate nodes at %s" % (len(_candidate_molecules), tp))
    
    tmp = changed_molecules.loc[changed_molecules.time_point == tp_next, :]
    _next_changed_molecules = sorted(list(set(tmp.uniprot_accession)))
    print("%d unique target nodes at %s" % (len(_next_changed_molecules), tp_next))
    for node in _candidate_molecules:
        if node in G:
            ds = sorted(downstream_nodes_direct(G, node))
            _changed = [x for x in ds if x in _next_changed_molecules]
            _other = sorted(list(set(ds) - set(_changed)))
        else:
            ds = []
            _changed = []
        if len(ds) > 0:
            ratio = len(_changed) / len(ds)
            other_nodes = pd.concat([
                other_nodes, pd.DataFrame({
                    'Node_ID': map(lambda x: x + "_" + tp_next, _other),
                    'uniprot_accession': _other,
                    'time_point': tp_next
                })
            ])
        else:
            ratio = -1

        r = pd.DataFrame([{'uniprot_accession': node, 'time_point': tp, 'Ratio': ratio,
                          'changed': len(_changed), 'total': len(ds),
                          'changed_nodes': ','.join(_changed), 'other_nodes': ','.join([x for x in ds if x not in _changed])
                          }])
        node_idx = "%s_%s" % (node, tp)
        impact = pd.concat([impact, r.rename(index = {0: node_idx})])
        #changed_idx = [x in candidate_proteins for x in candidates]
other_nodes = other_nodes.drop_duplicates().set_index('Node_ID')
other_nodes.to_csv(param.NORMAL_NODE_FILE, sep = '\t')

stat_results = pd.DataFrame(
    index = [],
    columns = ['OddsRatio', 'pValue', 'FDR_BH', 'pop_changed', 'pop_total']
    )

for ptr in range(len(param.TIME_POINTS[:-1])):
    tp = param.TIME_POINTS[ptr]
    next_tp = param.TIME_POINTS[ptr + 1]
    print("time point = %s" % tp)
    pop_total_nodes = len(molecules)
    pop_changed_molecules = changed_molecules.uniprot_accession.loc[changed_molecules.time_point == next_tp]
    pop_changed_molecules = list(set(pop_changed_molecules))
    pop_changed_nodes = len(pop_changed_molecules)
    pop_other_nodes = pop_total_nodes - pop_changed_nodes
    print("> population nodes = %d" % pop_total_nodes)
    print("> population changed nodes = %d" % pop_changed_nodes)
    _df = impact.loc[impact.time_point == tp, :].loc[impact.Ratio >= 0, :]
    _stats = pd.DataFrame(index = _df.index, columns = stat_results.columns)
    _stats.loc[:, 'pop_changed'] = pop_changed_nodes
    _stats.loc[:, 'pop_total'] = pop_total_nodes
    for idx in _df.index:
        pathway_changed_nodes = _df.changed[idx]
        pathway_other_nodes = _df.total[idx] - pathway_changed_nodes
        tab = [[pathway_changed_nodes, pathway_other_nodes],
                [pop_changed_nodes - pathway_changed_nodes, pop_other_nodes - pathway_other_nodes]]
        print(_df)
        odds_ratio, pvalue = stats.fisher_exact(tab)
        _stats.loc[idx, ['OddsRatio', 'pValue']] = pd.Series({'OddsRatio': odds_ratio, 'pValue': pvalue})
    stat_results = pd.concat([stat_results, _stats])

pval_corr = multi.multipletests(stat_results.pValue, alpha = 0.05, method = 'fdr_bh')
stat_results.loc[:, 'FDR_BH'] = pval_corr[1]
print(impact.shape)
print(stat_results.shape)
pd.concat([impact, stat_results], axis = 1).to_csv(param.IMPACT_RESULT_FILE, sep = '\t')
