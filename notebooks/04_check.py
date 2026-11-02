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

# 1. scaffold 다양성
def scaf(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else None
    except Exception:
        return None

uniq = g.drop_duplicates("Ligand InChI Key").copy()
uniq["scaffold"] = uniq["Ligand SMILES"].apply(scaf)
ok = uniq[uniq["scaffold"].notna()]

print("=" * 60)
print(f"G12D 단독 IC50(=) 측정: {len(g):,}")
print(f"고유 분자: {len(uniq):,}  (SMILES 파싱 성공 {len(ok):,})")
print(f"고유 scaffold: {ok['scaffold'].nunique():,}")
vc = ok["scaffold"].value_counts()
print(f"분자 1개뿐인 scaffold: {(vc == 1).sum():,}")
print(f"상위 10개 scaffold가 덮는 분자 비율: {vc.head(10).sum() / len(ok):.1%}")
print("\n상위 10 scaffold 크기:", vc.head(10).tolist())

# 2. pIC50 분포
print("\n" + "=" * 60)
print("[pIC50 분포]")
print(g["pIC50"].describe().to_string())
print("\n히스토그램:")
h, edges = np.histogram(g["pIC50"], bins=12)
for c, lo, hi in zip(h, edges[:-1], edges[1:]):
    print(f"  {lo:5.1f}~{hi:5.1f}  {'#' * int(c / max(h) * 50)} {c}")

# 3. 특허/문헌 집중도
print("\n" + "=" * 60)
print("[특허별 분자 수 상위 10]")
print(g.groupby("Patent Number")["Ligand InChI Key"].nunique()
       .sort_values(ascending=False).head(10).to_string()
      if "Patent Number" in g.columns else "(Patent Number 컬럼 없음)")

uniq.to_csv("data/processed/g12d_unique.csv", index=False)
g.to_csv("data/processed/g12d_measurements.csv", index=False)