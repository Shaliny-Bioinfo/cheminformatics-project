import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
import plotly.express as px

# === USER SET INPUT FILE NAME ===

INPUT_FILENAME = "golden_triangle_infantshaliny.csv"

# === PATH SETUP ===
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROBLEM_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROBLEM_DIR, "data")
RESULTS_DIR = os.path.join(PROBLEM_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Build full path for input file
INPUT_FILE = os.path.join(DATA_DIR, INPUT_FILENAME)
if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"❌ Input file not found: {INPUT_FILE}")
print(f"📂 Using input file: {INPUT_FILE}")

# === READ CSV WITH ENCODING FIX ===
df = pd.read_csv(INPUT_FILE, encoding='latin1')
if 'SMILES' not in df.columns:
    raise ValueError("❌ CSV must contain a 'SMILES' column!")

# === COMPUTE DESCRIPTORS ===
MW_values, LogP_values, TPSA_values, HBD_values, HBA_values, RO5_violations = [], [], [], [], [], []

for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        MW = Descriptors.MolWt(mol)
        LogP = Crippen.MolLogP(mol)
        TPSA = rdMolDescriptors.CalcTPSA(mol)
        HBD = Lipinski.NumHDonors(mol)
        HBA = Lipinski.NumHAcceptors(mol)

        # Rule of Five violations
        RO5_violation = 0
        if MW > 500: RO5_violation += 1
        if LogP > 5: RO5_violation += 1
        if HBD > 5: RO5_violation += 1
        if HBA > 10: RO5_violation += 1

        MW_values.append(MW)
        LogP_values.append(LogP)
        TPSA_values.append(TPSA)
        HBD_values.append(HBD)
        HBA_values.append(HBA)
        RO5_violations.append(RO5_violation)
    else:
        MW_values.append(None)
        LogP_values.append(None)
        TPSA_values.append(None)
        HBD_values.append(None)
        HBA_values.append(None)
        RO5_violations.append(None)

df['MW'] = MW_values
df['LogP'] = LogP_values
df['TPSA'] = TPSA_values
df['HBD'] = HBD_values
df['HBA'] = HBA_values
df['RO5_violations'] = RO5_violations

# === FILTER FLAGS ===
df['MW_Pass'] = df['MW'] <= 350
df['LogP_Pass'] = df['LogP'] <= 3.5
df['TPSA_Pass'] = df['TPSA'] <= 90
df['RO5_Pass'] = df['RO5_violations'] == 0
df['GoldenTriangle_Pass'] = df['MW_Pass'] & df['LogP_Pass'] & df['TPSA_Pass'] & df['RO5_Pass']

# === SUMMARY ===
total = len(df)
print(f"\nTotal compounds: {total}")
print(f"Pass MW filter: {df['MW_Pass'].sum()}")
print(f"Pass LogP filter: {df['LogP_Pass'].sum()}")
print(f"Pass TPSA filter: {df['TPSA_Pass'].sum()}")
print(f"Pass RO5 filter: {df['RO5_Pass'].sum()}")
print(f"Pass all 4 filters (Lead-like): {df['GoldenTriangle_Pass'].sum()}\n")

# === SAVE FILTERED CSV & HTML (with HBD, HBA, RO5 violations included) ===
base_name = os.path.splitext(INPUT_FILENAME)[0]
filtered_df = df[df['GoldenTriangle_Pass']].copy()
filtered_columns = ['SMILES','MW','LogP','TPSA','HBD','HBA','RO5_violations']
filtered_df = filtered_df[filtered_columns]

csv_output = os.path.join(RESULTS_DIR, f"{base_name}_Filtered.csv")
html_output = os.path.join(RESULTS_DIR, f"{base_name}_Filtered.html")

filtered_df.to_csv(csv_output, index=False)
filtered_df.to_html(html_output, index=False)

print(f"✅ Filtered CSV saved: {csv_output}")
print(f"✅ Filtered HTML table saved: {html_output}")

# === COLOR-CODE FOR PLOT ===
def assign_color(row):
    if row['GoldenTriangle_Pass']:
        return "Green: Lead-like"
    fail_flags = []
    if not row['MW_Pass']: fail_flags.append("MW")
    if not row['LogP_Pass']: fail_flags.append("LogP")
    if not row['TPSA_Pass']: fail_flags.append("TPSA")
    if not row['RO5_Pass']: fail_flags.append("RO5")
    if len(fail_flags) == 0:
        return "Grey: Unknown"
    return "Red: Fail " + "/".join(fail_flags)

df['ColorStatus'] = df.apply(assign_color, axis=1)

# === VISUALIZE WITH FULL HOVER INFO ===
hover_text = []
for idx, row in df.iterrows():
    hover_text.append(
        f"SMILES: {row['SMILES']}<br>"
        f"MW: {row['MW']:.2f}<br>"
        f"LogP: {row['LogP']:.2f}<br>"
        f"TPSA: {row['TPSA']:.2f}<br>"
        f"HBD: {row['HBD']}<br>"
        f"HBA: {row['HBA']}<br>"
        f"RO5 violations: {row['RO5_violations']}<br>"
        f"Lead-like: {'Yes' if row['GoldenTriangle_Pass'] else 'No'}"
    )

plot_output = os.path.join(RESULTS_DIR, f"{base_name}_Plot.html")
fig = px.scatter(
    df,
    x="MW",
    y="LogP",
    color="ColorStatus",
    title="Golden Triangle Lead-like Filtering",
    labels={"MW": "Molecular Weight", "LogP": "LogP"},
    hover_name=df['SMILES'],
    hover_data={'hover_text': hover_text},
    color_discrete_map={
        "Green: Lead-like": "green",
        "Red: Fail MW": "red",
        "Red: Fail LogP": "blue",
        "Red: Fail TPSA": "orange",
        "Red: Fail RO5": "purple",
        "Grey: Unknown": "grey"
    }
)
fig.update_traces(hovertemplate="%{customdata[0]}")
fig.write_html(plot_output)
print(f"📊 Interactive color-coded plot with full hover info saved: {plot_output}")
