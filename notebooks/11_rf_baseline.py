import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

SEEDS = [0, 1, 2, 3, 4]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def fingerprint(smi):
    m = Chem.MolFromSmiles(smi)
    return np.array(gen.GetFingerprint(m), dtype=np.uint8) if m else None

def scaffold_split(df, test_frac=0.2, seed=42):
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

def evaluate(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "MAE": mean_absolute_error(y_true, y_pred),
    }

rows = []

for name in ["cellular", "biochem", "mixed"]:
    df = pd.read_csv(f"data/processed/dataset_{name}.csv")
    fps = df["Ligand SMILES"].apply(fingerprint)
    df = df[fps.notna()].reset_index(drop=True)
    X = np.stack(fps[fps.notna()].values)
    y = df["pIC50"].values

    for seed in SEEDS:
        # scaffold split
        tr, te = scaffold_split(df, seed=seed)
        rf = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1)
        rf.fit(X[tr], y[tr])
        rows.append({"dataset": name, "split": "scaffold", "model": "RF",
                     "seed": seed, **evaluate(y[te], rf.predict(X[te]))})
        rows.append({"dataset": name, "split": "scaffold", "model": "mean",
                     "seed": seed,
                     **evaluate(y[te], np.full(len(te), y[tr].mean()))})

        # random split (비교용)
        rtr, rte = train_test_split(np.arange(len(y)), test_size=0.2,
                                    random_state=seed)
        rf2 = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1)
        rf2.fit(X[rtr], y[rtr])
        rows.append({"dataset": name, "split": "random", "model": "RF",
                     "seed": seed, **evaluate(y[rte], rf2.predict(X[rte]))})

res = pd.DataFrame(rows)
res.to_csv("data/processed/rf_baseline_raw.csv", index=False)

summary = (res.groupby(["dataset", "split", "model"])[["R2", "RMSE", "MAE"]]
             .agg(["mean", "std"]).round(3))
print(summary.to_string())