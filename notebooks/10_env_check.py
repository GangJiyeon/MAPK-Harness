import sys, platform
print(f"Python {sys.version.split()[0]} / {platform.machine()}")

mods = {}

try:
    import pandas as pd; mods["pandas"] = pd.__version__
except Exception as e: mods["pandas"] = f"FAIL {e}"

try:
    import numpy as np; mods["numpy"] = np.__version__
except Exception as e: mods["numpy"] = f"FAIL {e}"

try:
    import scipy; mods["scipy"] = scipy.__version__
    from scipy.integrate import solve_ivp
    r = solve_ivp(lambda t, y: -0.5 * y, [0, 10], [1.0], dense_output=True)
    mods["solve_ivp"] = f"ok, y(10)={r.y[0][-1]:.4f} (expect ~0.0067)"
except Exception as e: mods["scipy"] = f"FAIL {e}"

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    import rdkit
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mods["rdkit"] = f"{rdkit.__version__}, scaffold={MurckoScaffold.MurckoScaffoldSmiles(mol=m)}"
except Exception as e: mods["rdkit"] = f"FAIL {e}"

try:
    import sklearn
    from sklearn.ensemble import RandomForestRegressor
    mods["sklearn"] = sklearn.__version__
except Exception as e: mods["sklearn"] = f"FAIL {e}"

try:
    import torch
    mods["torch"] = f"{torch.__version__}, mps={torch.backends.mps.is_available()}"
except Exception as e: mods["torch"] = f"FAIL {e}"

try:
    import torch_geometric
    from torch_geometric.nn import GCNConv
    from torch_geometric.data import Data
    import torch
    d = Data(x=torch.randn(4, 3), edge_index=torch.tensor([[0,1,2],[1,2,3]]))
    out = GCNConv(3, 8)(d.x, d.edge_index)
    mods["torch_geometric"] = f"{torch_geometric.__version__}, GCNConv out={tuple(out.shape)}"
except Exception as e: mods["torch_geometric"] = f"FAIL {e}"

try:
    import matplotlib; mods["matplotlib"] = matplotlib.__version__
except Exception as e: mods["matplotlib"] = f"FAIL {e}"

try:
    import basico; mods["basico(COPASI)"] = "ok"
except Exception as e: mods["basico(COPASI)"] = f"FAIL {e}"

print("-" * 60)
for k, v in mods.items():
    print(f"{k:20} {v}")