import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import plotly.express as px

# === USER SET INPUT FILE NAME ===
# Just change this filename (must be in the 'data' folder)
INPUT_FILENAME = "chemical_space_infantshaliny.csv"

# === PATH SETUP ===
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROBLEM_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROBLEM_DIR, "data")
RESULTS_DIR = os.path.join(PROBLEM_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Build full input path
INPUT_FILE = os.path.join(DATA_DIR, INPUT_FILENAME)
if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"❌ Input file not found: {INPUT_FILE}")
print(f"📂 Using input file: {INPUT_FILE}")

# === READ CSV ===
df = pd.read_csv(INPUT_FILE, encoding='latin1')

# === COMPUTE DESCRIPTORS ===
MW_values, LogP_values, TPSA_values = [], [], []
for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        MW_values.append(Descriptors.MolWt(mol))
        LogP_values.append(Crippen.MolLogP(mol))
        TPSA_values.append(rdMolDescriptors.CalcTPSA(mol))
    else:
        MW_values.append(None)
        LogP_values.append(None)
        TPSA_values.append(None)

df['MW'] = MW_values
df['LogP'] = LogP_values
df['TPSA'] = TPSA_values

# === FLAG DRUG-LIKE / OUTLIERS ===
df['Type'] = 'Other'
df.loc[(df['MW'].between(200, 350)) & (df['LogP'].between(0, 3.5)), 'Type'] = 'Drug-like'

# Select 4 largest MW as outliers for example
outliers_idx = df['MW'].nlargest(4).index
df.loc[outliers_idx, 'Type'] = 'Outlier'

# === SAVE CSV ===
csv_output = os.path.join(RESULTS_DIR, f"{os.path.splitext(INPUT_FILENAME)[0]}_Descriptors.csv")
df.to_csv(csv_output, index=False)

# === INTERACTIVE SCATTER PLOT WITH HOVER INFO ===
hover_text = [
    f"SMILES: {row['SMILES']}<br>"
    f"MW: {row['MW']:.2f}<br>"
    f"LogP: {row['LogP']:.2f}<br>"
    f"TPSA: {row['TPSA']:.2f}<br>"
    f"Type: {row['Type']}"
    for idx, row in df.iterrows()
]

fig = px.scatter(
    df,
    x="MW",
    y="LogP",
    color="Type",
    hover_data={'hover_text': hover_text},
    hover_name=df['SMILES'],
    title="Chemical Space: MW vs LogP",
    color_discrete_map={"Drug-like":"green", "Outlier":"red", "Other":"blue"}
)
fig.update_traces(hovertemplate="%{customdata[0]}")

plot_output = os.path.join(RESULTS_DIR, f"{os.path.splitext(INPUT_FILENAME)[0]}_Plot.html")
fig.write_html(plot_output)

print(f"✅ Problem 2 outputs saved in {RESULTS_DIR}")
