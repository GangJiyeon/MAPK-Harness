import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

df = pd.read_csv("data/processed/g12d_assay_classified.csv", dtype=str)
df["pIC50"] = pd.to_numeric(df["pIC50"], errors="coerce")

# 제외
df = df[~df["assay_class"].isin(["biochem_unspecified", "unknown"])].copy()
print(f"제외 후: {len(df):,}")

# 같은 분자가 같은 기전에서 여러 번 측정됐으면 중위값
def dedupe(sub):
    return (sub.groupby(["Ligand InChI Key", "Ligand SMILES", "Patent Number"],
                        as_index=False)
               .agg(pIC50=("pIC50", "median"), n=("pIC50", "size")))

sets = {
    "cellular": df[df["assay_class"] == "cellular"],
    "biochem":  df[df["assay_class"].isin(["exchange", "binding"])],
    "mixed":    df,
}

def scaf(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else None
    except Exception:
        return None

for name, sub in sets.items():
    d = dedupe(sub)
    d["scaffold"] = d["Ligand SMILES"].apply(scaf)
    d = d[d["scaffold"].notna()]
    print(f"\n[{name}] 분자 {len(d):,} / scaffold {d['scaffold'].nunique():,} "
          f"/ 특허 {d['Patent Number'].nunique()}")
    print(f"  pIC50 평균 {d['pIC50'].mean():.2f} 표준편차 {d['pIC50'].std():.2f}")
    d.to_csv(f"data/processed/dataset_{name}.csv", index=False)