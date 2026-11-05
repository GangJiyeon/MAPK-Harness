"""데이터 정제 과정의 확인 스크립트.

1. pIC50 상/하한 이상치가 실제로 어떤 값인지
2. biochem_unspecified 422건이 어떤 assay를 측정한 것인지
"""
import pandas as pd

# --- 1. pIC50 이상치 확인 ---
g = pd.read_csv("data/processed/g12d_measurements.csv")

print("[상한 이상치: pIC50 > 10.6]")
print(g[g["pIC50"] > 10.6][
    ["BindingDB Ligand Name", "ic50_nM", "pIC50", "Curation/DataSource"]
].to_string())
# → 2개 특허에 집중, 연속 화합물인데 값이 불규칙 → 측정 한계 초과로 판단, 제외

print("\n[하한 확인: pIC50 < 4.4]")
g = g[g["pIC50"] <= 10].copy()
print(g[g["pIC50"] < 4.4][
    ["BindingDB Ligand Name", "ic50_nM", "pIC50"]
].to_string())
# → 50~190 µM은 측정 가능 범위. 약한 화합물도 학습에 필요하므로 유지

# --- 2. biochem_unspecified가 무엇을 측정했는지 ---
df = pd.read_csv("data/processed/g12d_assay_classified.csv", dtype=str)
u = df[df["assay_class"] == "biochem_unspecified"]

print(f"\n[biochem_unspecified] {len(u)}건, "
      f"고유 DESCRIPTION {u['DESCRIPTION'].nunique()}개")
for d in u["DESCRIPTION"].dropna().unique()[:1]:
    print(str(d)[:600])
# → 전부 동일한 assay. covalent labeling(공유결합 표지)이며
#   본문에 K-Ras(1-169) G12C로 명시 → G12D 라벨과 불일치. 제외