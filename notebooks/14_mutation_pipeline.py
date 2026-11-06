import pandas as pd
import re

df = pd.read_csv("data/processed/kras_classified.csv", dtype=str)
df["ic50_nM"] = pd.to_numeric(df["ic50_nM"], errors="coerce")

# assay 설명 조인
rsid = pd.read_csv("data/raw/BindingDB_rsid_eaids.tsv", sep="\t", dtype=str)
rsid = rsid.rename(columns={"REACTANT_SET_ID": "BindingDB Reactant_set_id"})
assays = pd.read_csv("data/raw/BindingDB_Assays.tsv", sep="\t", dtype=str,
                     on_bad_lines="skip")

m = df.merge(rsid, on="BindingDB Reactant_set_id", how="left")
m[["ENTRYID", "ASSAYID"]] = m["ENTRYID_ASSAYID"].str.split("_", n=1, expand=True)
m = m.merge(assays, on=["ENTRYID", "ASSAYID"], how="left")

def classify(s):
    if any(k in s for k in ["phosphorylation of erk", "perk", "cell",
                            "proliferation", "viability"]):
        return "cellular"
    if any(k in s for k in ["covalent", "covalently label"]):
        return "covalent"
    if any(k in s for k in ["nucleotide exchange", "gdp exchange",
                            "sos1-mediated", "sos-mediated", "bodipy"]):
        return "exchange"
    if any(k in s for k in ["binding", "tr-fret", "trfret", "htrf"]):
        return "binding"
    if "biochemical" in s:
        return "biochem_unspecified"
    return "unknown"

both = (m["ASSAY_NAME"].fillna("") + " | " + m["DESCRIPTION"].fillna("")).str.lower()
m["assay_class"] = both.apply(classify)

# 파이프라인 통과분만
ok = m[(m["is_complex"] == "False") & (m["ic50_rel"] == "=") &
       m["ic50_nM"].notna() & (m["ic50_nM"] > 0)].copy()
ok["pIC50"] = 9 - pd.np.log10(ok["ic50_nM"]) if False else (9 - __import__("numpy").log10(ok["ic50_nM"]))
ok = ok[ok["pIC50"] <= 10]

muts = ["G12D", "G12C", "G12V", "WT", "G13D"]

print("[변이 × assay 기전, 고유분자 수]")
t = (ok[ok["mutation"].isin(muts)]
     .groupby(["mutation", "assay_class"])["Ligand InChI Key"]
     .nunique().unstack(fill_value=0))
print(t.to_string())

# 공통 분자 — 선택성 계산 가능 여부
print("\n[변이 간 공통 분자 수 (전체 기전)]")
sets = {mu: set(ok[ok["mutation"] == mu]["Ligand InChI Key"]) for mu in muts}
for a in muts:
    row = [f"{len(sets[a] & sets[b]):5,}" for b in muts]
    print(f"  {a:6} " + " ".join(row))
print("  " + " " * 6 + " ".join(f"{mu:>5}" for mu in muts))

print("\n[cellular 기전만, 공통 분자]")
cel = ok[ok["assay_class"] == "cellular"]
csets = {mu: set(cel[cel["mutation"] == mu]["Ligand InChI Key"]) for mu in muts}
print(f"  G12D ∩ WT   : {len(csets['G12D'] & csets['WT']):,}")
print(f"  G12D ∩ G12C : {len(csets['G12D'] & csets['G12C']):,}")
print(f"  G12D ∩ G12V : {len(csets['G12D'] & csets['G12V']):,}")

bio = ok[ok["assay_class"].isin(["exchange", "binding"])]
bsets = {mu: set(bio[bio["mutation"] == mu]["Ligand InChI Key"]) for mu in muts}
print(f"biochem G12D ∩ WT   : {len(bsets['G12D'] & bsets['WT']):,}")
print(f"biochem G12D ∩ G12C : {len(bsets['G12D'] & bsets['G12C']):,}")

# 같은 기전끼리만
for cls in ["exchange", "binding"]:
    s = ok[ok["assay_class"] == cls]
    a = set(s[s["mutation"] == "G12D"]["Ligand InChI Key"])
    b = set(s[s["mutation"] == "WT"]["Ligand InChI Key"])
    print(f"{cls:10} G12D ∩ WT: {len(a & b):,}")

ex = ok[ok["assay_class"] == "exchange"]
g = ex[ex["mutation"] == "G12D"].groupby("Ligand InChI Key")["pIC50"].median()
w = ex[ex["mutation"] == "WT"].groupby("Ligand InChI Key")["pIC50"].median()
common = g.index.intersection(w.index)
sel = (g[common] - w[common]).sort_values(ascending=False)

print(f"공통 {len(sel)}개")
print(sel.describe().to_string())
print("\n[선택성 상위 5]")
print(sel.head().to_string())
print("\n[선택성 하위 5]")
print(sel.tail().to_string())

import numpy as np
h, e = np.histogram(sel, bins=10)
print("\n[분포]")
for c, lo, hi in zip(h, e[:-1], e[1:]):
    print(f"  {lo:5.1f}~{hi:5.1f}  {'#' * int(c/max(h)*40)} {c}")


# 선택성과 G12D pIC50을 한 표로
tbl = pd.DataFrame({
    "pIC50_G12D": g[common],
    "pIC50_WT": w[common],
    "selectivity": sel,
})

print(f"pIC50_G12D와 selectivity 상관계수: "
      f"{tbl['pIC50_G12D'].corr(tbl['selectivity']):.3f}")

top_potency = set(tbl.nlargest(5, "pIC50_G12D").index)
top_select  = set(tbl.nlargest(5, "selectivity").index)
print(f"\npIC50 Top5 ∩ 선택성 Top5: {len(top_potency & top_select)}개")

print("\n[pIC50 Top 10]")
print(tbl.nlargest(10, "pIC50_G12D").round(2).to_string())

print("\n[선택성 Top 10]")
print(tbl.nlargest(10, "selectivity").round(2).to_string())

tbl = tbl.reset_index()
smi = ok.drop_duplicates("Ligand InChI Key").set_index("Ligand InChI Key")["Ligand SMILES"]
tbl["SMILES"] = tbl["Ligand InChI Key"].map(smi)
tbl.to_csv("data/processed/selectivity_g12d_wt.csv", index=False)
print(tbl.shape, tbl["SMILES"].notna().sum())