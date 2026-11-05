import pandas as pd
df = pd.read_csv("data/processed/g12d_assay_classified.csv", dtype=str)
u = df[df["assay_class"] == "biochem_unspecified"]
for d in u["DESCRIPTION"].dropna().unique()[:3]:
    print("-" * 60)
    print(str(d)[:600])

    