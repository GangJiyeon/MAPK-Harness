import pandas as pd

g = pd.read_csv("data/processed/g12d_measurements.csv", dtype=str)
print(f"G12D 측정 {len(g):,}건")

rsid = pd.read_csv("data/raw/BindingDB_rsid_eaids.tsv", sep="\t", dtype=str)
print(f"\nrsid 파일 컬럼: {list(rsid.columns)}")
print(rsid.head(3).to_string())

assays = pd.read_csv("data/raw/BindingDB_Assays.tsv", sep="\t", dtype=str,
                     on_bad_lines="skip")
print(f"\nassays 파일 컬럼: {list(assays.columns)}")
print(assays.head(3).to_string())



# 조인
rsid = rsid.rename(columns={"REACTANT_SET_ID": "BindingDB Reactant_set_id"})
m = g.merge(rsid, on="BindingDB Reactant_set_id", how="left")
print(f"\nrsid 매칭: {m['ENTRYID_ASSAYID'].notna().sum():,} / {len(m):,}")

m[["ENTRYID", "ASSAYID"]] = m["ENTRYID_ASSAYID"].str.split("_", n=1, expand=True)
m2 = m.merge(assays, on=["ENTRYID", "ASSAYID"], how="left")
print(f"assay 설명 매칭: {m2['DESCRIPTION'].notna().sum():,} / {len(m2):,}")

if m2["DESCRIPTION"].notna().sum() > 0:
    print("\n[ASSAY_NAME 상위 20]")
    print(m2["ASSAY_NAME"].value_counts().head(20).to_string())
    print("\n[DESCRIPTION 샘플 3개]")
    for d in m2["DESCRIPTION"].dropna().unique()[:3]:
        print("-" * 50)
        print(str(d)[:400])
    m2.to_csv("data/processed/g12d_with_assay.csv", index=False)
else:
    print("\n특허 데이터에 assay 설명 없음 → 한계로 기록")