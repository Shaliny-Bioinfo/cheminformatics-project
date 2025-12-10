import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import plotly.graph_objects as go
import os

# === Input and Output paths ===
input_csv = "../data/fda_drugs_infantshaliny.csv"  # <-- change if your file name differs
output_html = "../results/fda_drugs_farrin_RO5_Violations_Interactive.html"
output_csv = "../results/fda_drugs_farrin_RO5_Analysis.csv"

# === Read CSV ===
df = pd.read_csv(input_csv)

# Check available columns
print("Columns in input CSV:", df.columns.tolist())

# Try to locate the SMILES column automatically
smiles_col = None
for col in df.columns:
    if 'smiles' in col.lower():
        smiles_col = col
        break

if smiles_col is None:
    raise KeyError("No column containing SMILES found in the CSV file.")

# === Calculate Lipinski properties ===
mw_list, logp_list, hdonor_list, hacceptor_list = [], [], [], []

for smi in df[smiles_col]:
    if isinstance(smi, str) and smi.strip():
        mol = Chem.MolFromSmiles(smi)
    else:
        mol = None

    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hdonor = Descriptors.NumHDonors(mol)
        hacceptor = Descriptors.NumHAcceptors(mol)
    else:
        mw, logp, hdonor, hacceptor = [None]*4

    mw_list.append(mw)
    logp_list.append(logp)
    hdonor_list.append(hdonor)
    hacceptor_list.append(hacceptor)

df["MW"] = mw_list
df["LogP"] = logp_list
df["HDonors"] = hdonor_list
df["HAcceptors"] = hacceptor_list

# === Lipinski violations ===
df["RO5_violations"] = (
    (df["MW"] > 500).astype(int)
    + (df["LogP"] > 5).astype(int)
    + (df["HDonors"] > 5).astype(int)
    + (df["HAcceptors"] > 10).astype(int)
)

# === Count number of compounds violating each rule ===
violation_counts = {
    "Molecular Weight (>500)": (df["MW"] > 500).sum(),
    "LogP (>5)": (df["LogP"] > 5).sum(),
    "H Donors (>5)": (df["HDonors"] > 5).sum(),
    "H Acceptors (>10)": (df["HAcceptors"] > 10).sum()
}

# === Create bar chart ===
fig = go.Figure(
    data=[go.Bar(
        x=list(violation_counts.keys()),
        y=list(violation_counts.values()),
        text=[f"{v} compounds" for v in violation_counts.values()],
        textposition='auto',
        marker_color=['#FF6F61', '#6B5B95', '#88B04B', '#F7CAC9']
    )]
)

fig.update_layout(
    title="Lipinski Rule of Five Violations (FDA Drugs)",
    xaxis_title="Property",
    yaxis_title="Number of Compounds Violating",
    template="plotly_dark"
)

# === Save outputs ===
os.makedirs(os.path.dirname(output_html), exist_ok=True)
fig.write_html(output_html)
df.to_csv(output_csv, index=False)

print(f"✅ Interactive plot saved to: {output_html}")
print(f"📊 CSV analysis saved to: {output_csv}")
