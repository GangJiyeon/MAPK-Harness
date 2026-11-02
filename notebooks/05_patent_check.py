import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

df = pd.read_csv("data/processed/kras_classified.csv", dtype=str)
df["ic50_nM"] = pd.to_numeric(df["ic50_nM"], errors="coerce")

g = df[(df["mutation"] == "G12D") &
       (df["is_complex"] == "False") &
       (df["ic50_rel"] == "=") &
       df["ic50_nM"].notna() &
       (df["ic50_nM"] > 0)].copy()
g["pIC50"] = 9 - np.log10(g["ic50_nM"])
g = g[g["pIC50"] <= 10].copy()

uniq = g.drop_duplicates("Ligand InChI Key").copy()

def scaf(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else None
    except Exception:
        return None

uniq["scaffold"] = uniq["Ligand SMILES"].apply(scaf)
uniq = uniq[uniq["scaffold"].notna()]

print("=" * 60)
print(f"최종 측정 {len(g):,}건 / 고유 분자 {len(uniq):,} / scaffold {uniq['scaffold'].nunique():,}")

print("\n[특허 수]")
print(f"고유 특허: {uniq['Patent Number'].nunique():,}")

vc = uniq["Patent Number"].value_counts()
print(f"\n상위 10 특허가 덮는 분자 비율: {vc.head(10).sum()/len(uniq):.1%}")
print("\n[특허별 분자 수 상위 15]")
print(vc.head(15).to_string())

print("\n[특허별 분자 수 분포]")
print(vc.describe().to_string())

# 특허 하나가 scaffold 하나인지
print("\n[특허별 고유 scaffold 수 상위 15]")
ps = uniq.groupby("Patent Number")["scaffold"].nunique().sort_values(ascending=False)
print(ps.head(15).to_string())
print(f"\nscaffold 1개뿐인 특허: {(ps == 1).sum():,} / {len(ps):,}")

# scaffold가 여러 특허에 걸쳐 있는지 (겹침 = split 누수 위험)
sp = uniq.groupby("scaffold")["Patent Number"].nunique()
print(f"\n2개 이상 특허에 등장하는 scaffold: {(sp > 1).sum():,} / {len(sp):,}")

uniq.to_csv("data/processed/g12d_unique.csv", index=False)
g.to_csv("data/processed/g12d_measurements.csv", index=False)