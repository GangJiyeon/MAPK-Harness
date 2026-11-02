import pandas as pd
import re

df = pd.read_csv("data/processed/kras_raw.csv", dtype=str)

name = df["Target Name"].fillna("")
low = name.str.lower()

# 변이 추출
def mutation(s):
    muts = re.findall(r'\b([GQ]\d{2}[A-Z])\b', s.upper())
    muts = [m for m in muts if m not in ("C51S", "C80L", "C118S", "C118A")]
    return ",".join(sorted(set(muts))) if muts else "WT"

df["mutation"] = name.apply(mutation)

# 복합체 여부 (파트너 단백질과 함께 측정)
partners = ["sevenless", "sos", "raf", "kinase", "isomerase", "pi3k",
            "phosphatidylinositol"]
df["is_complex"] = low.apply(lambda s: any(p in s for p in partners))

# 값 파싱: 관계기호 분리
def parse(v):
    if pd.isna(v):
        return pd.Series([None, None])
    v = str(v).strip()
    m = re.match(r'^([<>=~]*)\s*([\d.eE+-]+)$', v)
    if not m:
        return pd.Series([None, None])
    rel = m.group(1) or "="
    try:
        return pd.Series([rel, float(m.group(2))])
    except ValueError:
        return pd.Series([None, None])

df[["ic50_rel", "ic50_nM"]] = df["IC50 (nM)"].apply(parse)

print("=" * 60)
print("[변이별]")
print(df["mutation"].value_counts().head(15).to_string())

print("\n[복합체 여부]")
print(df["is_complex"].value_counts().to_string())

print("\n[IC50 관계기호]")
print(df["ic50_rel"].value_counts(dropna=False).to_string())

g12d = df[(df["mutation"] == "G12D")]
print("\n" + "=" * 60)
print(f"G12D 전체: {len(g12d):,}")
for flag in [False, True]:
    sub = g12d[g12d["is_complex"] == flag]
    label = "복합체" if flag else "단독"
    ic = sub[(sub["ic50_rel"] == "=") & sub["ic50_nM"].notna()]
    print(f"  {label}: {len(sub):,}행, IC50(=) {len(ic):,}건, "
          f"고유분자 {ic['Ligand InChI Key'].nunique():,}")

print("\n[G12D 단독, 출처별]")
solo = g12d[~g12d["is_complex"]]
print(solo["Curation/DataSource"].value_counts().to_string())

print("\n[G12D 단독, Target Name 표기]")
print(solo["Target Name"].value_counts().head(20).to_string())

df.to_csv("data/processed/kras_classified.csv", index=False)