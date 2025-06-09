import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor

def get_beads(cas):
    """Return list of beads for a given CAS, trying .text/.txt under Working & Half_working."""
    for root in ["Working","Half_working"]:
        for ext in [".text",".txt"]:
            p = os.path.join(root, cas, f"{cas}{ext}")
            if os.path.isfile(p):
                with open(p) as f:
                    lines = f.readlines()[6:]
                beads = []
                for L in lines:
                    if not L.strip(): break
                    parts = L.split()
                    if len(parts) > 1:
                        beads.append(parts[1])
                return beads
    return []

# 1) Load errors
df = pd.read_csv("delta_testing_653.csv")
df["abs_error"] = (df.iloc[:,1] - df.iloc[:,2]).abs()

# 2) Build bead vocabulary and per‐CAS bead lists
vocab = set()
cas_beads = {}
for cas in df.iloc[:,0]:
    b = get_beads(str(cas))
    cas_beads[cas] = b
    vocab.update(b)
vocab = sorted(vocab)

# 3) Create feature matrix X and target y
X = np.array([[cas_beads[cas].count(bead) for bead in vocab] for cas in df.iloc[:,0]])
y = df["abs_error"].values

# 4) Fit Random Forest
rf = RandomForestRegressor(n_estimators=100, random_state=0)
rf.fit(X, y)

# 5) Extract importances
imp_df = pd.DataFrame({
    "bead": vocab,
    "importance": rf.feature_importances_
}).sort_values("importance", ascending=False)

# 6) Save and inspect top beads
imp_df.to_csv("bead_importance.csv", index=False)
print(imp_df.head(20))
