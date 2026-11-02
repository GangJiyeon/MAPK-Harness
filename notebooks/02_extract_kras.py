import pandas as pd

SRC = "data/raw/BindingDB_All.tsv"
OUT = "data/processed/kras_raw.csv"

COLS = [
    "BindingDB Reactant_set_id",    # assay 설명 조인키
    "Ligand SMILES",                # 모델입력
    "Ligand InChI Key",             # 중복판정
    "BindingDB Ligand Name",        
    "Target Name",                  # G12D/G12C 들어 있을 것??
    "Target Source Organism According to Curator or DataSource",    # human 필터링
    "Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)",                 # 측정값(mainly IC50)
    "Curation/DataSource",          # ChEMBL 유래인지 구분 > 중복제거
    "Article DOI", "PMID",          # 문헌 ID, 중복판정
    "Patent Number",
    "ChEMBL ID of Ligand",          # ChEMBL 대조용
    "UniProt (SwissProt) Primary ID of Target Chain 1",             # KRAS == PO1116
]

KEYS = ["kras", "k-ras", "gtpase kras", "p21 ras"]

chunks = []
total = 0
reader = pd.read_csv(
    SRC, sep="\t", chunksize=200_000, usecols=COLS,
    dtype=str, on_bad_lines="skip", engine="c", quoting=3,
)
for i, ch in enumerate(reader):
    total += len(ch)
    name = ch["Target Name"].fillna("").str.lower()
    mask = False
    for k in KEYS:
        mask = mask | name.str.contains(k, regex=False)
    hit = ch[mask]
    if len(hit):
        chunks.append(hit)
    if i % 10 == 0:
        print(f"chunk {i}, read {total:,}, hits {sum(len(c) for c in chunks):,}")

df = pd.concat(chunks, ignore_index=True)
df.to_csv(OUT, index=False)

print("=" * 50)
print(f"전체 행: {total:,}")
print(f"KRAS 행: {len(df):,}")
print("\n[Target Name 종류]")
print(df["Target Name"].value_counts().head(40).to_string())
print("\n[측정값 있는 행]")
for c in ["Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)"]:
    print(f"  {c}: {df[c].notna().sum():,}")
print("\n[데이터 출처]")
print(df["Curation/DataSource"].value_counts().to_string())