import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

RADIUS, NBITS = 2, 2048
gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)

def fingerprint(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return np.array(gen.GetFingerprint(m), dtype=np.uint8)

def scaffold_split(df, test_frac=0.2, seed=42):
    """골격 그룹을 무작위 순서로 나눔."""
    rng = np.random.default_rng(seed)
    groups = list(df.groupby("scaffold").indices.values())
    rng.shuffle(groups)
    n_train = int(len(df) * (1 - test_frac))
    train_idx, test_idx = [], []
    for idx in groups:
        if len(train_idx) + len(idx) <= n_train:
            train_idx.extend(idx)
        else:
            test_idx.extend(idx)
    return np.array(train_idx), np.array(test_idx)


for name in ["cellular", "biochem", "mixed"]:
    df = pd.read_csv(f"data/processed/dataset_{name}.csv")

    fps = df["Ligand SMILES"].apply(fingerprint)
    ok = fps.notna()
    df = df[ok].reset_index(drop=True)
    X = np.stack(fps[ok].values)
    y = df["pIC50"].values

    tr, te = scaffold_split(df)

    print(f"\n[{name}]")
    print(f"  분자 {len(df):,} / fingerprint {X.shape}")
    print(f"  학습 {len(tr):,} / 평가 {len(te):,}")
    print(f"  학습 scaffold {df.iloc[tr]['scaffold'].nunique()} / "
          f"평가 scaffold {df.iloc[te]['scaffold'].nunique()}")
    overlap = set(df.iloc[tr]['scaffold']) & set(df.iloc[te]['scaffold'])
    print(f"  겹치는 scaffold: {len(overlap)}  (0이어야 정상)")
    print(f"  평가셋 골격 크기: {df.iloc[te].groupby('scaffold').size().describe().to_dict()}")

    np.savez(f"data/processed/feat_{name}.npz",
             X=X, y=y, train_idx=tr, test_idx=te)


    