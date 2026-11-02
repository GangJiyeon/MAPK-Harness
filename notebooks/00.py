import pandas as pd
g = pd.read_csv("data/processed/g12d_measurements.csv")
print(g[g["pIC50"] > 10.6][["BindingDB Ligand Name", "ic50_nM", "pIC50", "Curation/DataSource"]].to_string())

g = g[g["pIC50"] <= 10].copy()
print(g[g["pIC50"] < 4.4][["BindingDB Ligand Name", "ic50_nM", "pIC50"]].to_string())