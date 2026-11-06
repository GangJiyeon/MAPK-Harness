import pandas as pd

df = pd.read_csv("data/processed/g12d_assay_classified.csv", dtype=str)
# ↑ 이건 G12D만 있으니 03 단계 파일로 다시 가야 함

df = pd.read_csv("data/processed/kras_classified.csv", dtype=str)
df["ic50_nM"] = pd.to_numeric(df["ic50_nM"], errors="coerce")

targets = ["G12D", "G12C", "G12V", "WT", "G13D"]

print("[변이별 파이프라인 통과 건수]")
for mut in targets:
    sub = df[df["mutation"] == mut]
    a = len(sub)
    b = len(sub[sub["is_complex"] == "False"])
    c = sub[(sub["is_complex"] == "False") & (sub["ic50_rel"] == "=")]
    c = c[c["ic50_nM"].notna() & (c["ic50_nM"] > 0)]
    print(f"  {mut:6} 전체 {a:6,} → 단독 {b:6,} → IC50(=) {len(c):6,} "
          f"→ 고유분자 {c['Ligand InChI Key'].nunique():5,}")