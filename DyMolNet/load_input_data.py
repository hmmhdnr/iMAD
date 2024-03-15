import pandas as pd
import os
import sys
import math

#
output_columns = [
    'Gene_ID', 'Group', 'UniProt_accession', 'Gene_Symbol', 'Gene_name',
    'time_point', 'ratio', 'pValue'
]

input_columns = {
    "group": "Group",
    "time_point": "time_point",
    "molecule_id": "Gene_ID",
    "gene_symbol": "Gene_symbol",
    "gene_name": "Gene_name",
    "value": "ratio",
    "p_value": "p_value",
    "q_value": "q_value"
}
#

def column_names(
        molecule_id = "Gene_ID",
        gene_symbol = "Gene_Symbol",
        gene_name = "Gene_name",
        value = "ratio",
        p_value = "p_value",
        q_value = "q_value"
):
    _dict = {
        "molecule_id": molecule_id,
        "gene_symbol": gene_symbol,
        "gene_name": gene_name,
        "value": value,
        "p_value": p_value,
        "q_value": q_value
    }
    return _dict

def _read_datasheet(path):
    fname, ext = os.path.splitext(path)
    if ext == ".tsv" or ext == ".txt":
        return pd.read_table(path, sep = "\t")
    elif ext == ".xlsx":
        return pd.read_excel(path)
    else:
        print("Error: datasheet was empty: %s" % path)
        sys.exit()

def input_proteome(input_columns, input_file = ""):
    _result = pd.DataFrame(
        index = [],
        columns = output_columns
    )

    if not os.path.exists(input_file):
        print("Error: couldn't read input datasheet %s" % input_file)
        sys.exit()
    _df = _read_datasheet(input_file)
    

    return _df

def input_RNAseq(
        input_file = "datasheet.txt",
        time_point = "t0",
        group_name = "disease",
        project_name = "",
        id_col = "Gene_ID",
        symbol_col = "symbol",
        uniprot_col = 'uniprot_accession',
        names_col = "Gene_name",
        values_col = "ratio",
        pvalues_col = "p_value",
        **kwargs        
):
    _result = pd.DataFrame(index = [], columns = [])

    if not os.path.exists(input_file):
        print("Error: couldn't read input datasheet %s" % input_file)
        return pd.DataFrame()
    _df = _read_datasheet(input_file)

    if id_col not in _df.columns:
        print("please specify valid column name for Molecule_ID: %s" % id_col)
        return pd.DataFrame()
    if values_col not in _df.columns:
        print("please specify valid column name for Expression_Values: %s" % values_col)
        return pd.DataFrame()
    if pvalues_col not in _df.columns:
        print("please specify valid column name for p_values: %s" % pvalues_col)
        return pd.DataFrame()

    _result["molecule_id"] = _df[id_col]
    _result["value"] = _df[values_col]
    _result["p_value"] = _df[pvalues_col]

    if symbol_col in _df.columns:
        _result["gene_symbol"] = _df[symbol_col]

    if names_col in _df.columns:
         _result["gene_names"] = _df[names_col]
    
    if uniprot_col in _df.columns:
        _result['uniprot_accession'] = _df[uniprot_col]

    return _result

## --> 2023/6/2 -->
def input_data():
    return pd.DataFrame()

def load_data(
        project_name, group_name, data_type, time_points,
        data_dir, id_col, uniprot_col, symbol_col,
        values_col, pvalues_col, **kwargs
    ):
    result = pd.DataFrame()
    if data_type == "mRNA_expression":
        _input_func = input_RNAseq
    elif data_type == "phospho_peptide":
        _input_func = input_proteome
    else:
        _input_func = input_data
    
    for time_point in time_points:
        _fname = "%s_%s_%s.tsv" % (project_name, group_name, time_point)
        _args = {
            "input_file": os.path.join(data_dir, _fname),
            "time_point": time_point,
            "group_name": group_name,
            "project_name": project_name,
            "id_col": id_col,
            "uniprot_col": uniprot_col,
            "symbol_col": symbol_col,
            "values_col": values_col,
            "pvalues_col": pvalues_col,
        }
        if "names_col" in kwargs.keys():
            _args["names_col"] = kwargs["names_col"]
        if "qvalues_col" in kwargs.keys():
            _args["qvalues_col"] = kwargs["qvalues_col"]
        if "case_col" in kwargs.keys():
            _args["case_col"] = kwargs["case_col"]
        if "control_col" in kwargs.keys():
            _args["control_col"] = kwargs["control_col"]

        _df = _input_func(**_args)
        if not _df.empty:
            _df = _df.set_index("molecule_id")
            _df["time_point"] = time_point
            result = pd.concat([result, _df])
    return result

# DEG selection by log2 expression ratio
def get_molecule_list_lfc(
        df,
        lfc_col = "ratio",
        lfc_threshold = 1.0,
        pvalues_col = "pValue",
        all_molecule_path = "all_molecules.txt",
        changed_molecule_path = "changed_molecules.txt"
    ):
    # choose candidate molecules
    molecules = [m for m in list(set(df.uniprot_accession)) if not pd.isnull(m)]
    changed = df.loc[df.loc[:, lfc_col] > lfc_threshold, :]
    changed = changed.loc[:, ["uniprot_accession", "time_point", "gene_symbol", "gene_names", lfc_col, pvalues_col]]
    changed["Node_ID"] = changed["uniprot_accession"].str.cat(changed["time_point"], sep = "_")
    changed_molecules = changed.drop_duplicates(subset = "Node_ID")
    changed_molecules = changed_molecules.dropna(subset = ["uniprot_accession"])
    changed_molecules = changed_molecules.set_index("Node_ID")
    ## output
    with open(all_molecule_path, "w") as f:
        f.write("\n".join(molecules))
    changed_molecules.to_csv(changed_molecule_path, sep = '\t')
    
    return [molecules, changed, changed_molecules]

# DEG selection by p-value
def get_molecule_list(
        df,
        pvalues_col = "pValue",
        p_threshold = 0.05,
        all_molecule_path = "all_molecules.txt",
        changed_molecule_path = "changed_molecules.txt"
    ):
    # choose candidate molecules
    molecules = [m for m in list(set(df.uniprot_accession)) if not pd.isnull(m)]
    changed = df.loc[df.loc[:, pvalues_col] < p_threshold, :]
    changed = changed.loc[:, ["uniprot_accession", "time_point", "gene_symbol", "gene_names", "value", pvalues_col]]
    changed["Node_ID"] = changed["uniprot_accession"].str.cat(changed["time_point"], sep = "_")
    changed_molecules = changed.drop_duplicates(subset = "Node_ID")
    changed_molecules = changed_molecules.dropna(subset = ["uniprot_accession"])
    changed_molecules = changed_molecules.set_index("Node_ID")
    ## output
    with open(all_molecule_path, "w") as f:
        f.write("\n".join(molecules))
    changed_molecules.to_csv(changed_molecule_path, sep = '\t')
    
    return [molecules, changed, changed_molecules]

if __name__ == "__main__":
    print("DyMolNet.load_input_data")
