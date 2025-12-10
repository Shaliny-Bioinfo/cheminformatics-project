#  Computational Drug Discovery using Cheminformatics Tools  
*A small-scale virtual screening and molecular property analysis pipeline built by **Infant Shaliny Arokia Lourdu Elango***  

---

##  Overview

This repository contains a **cheminformatics workflow** implemented in Python and RDKit to:

- Calculate basic molecular descriptors  
- Apply **Golden Triangle** virtual screening filters  
- Visualize **chemical space** using logP and TPSA  
- Assess **Lipinski’s Rule of Five (RO5)** compliance  

The workflow is demonstrated on a curated dataset of **20 small molecules** retrieved from **PubChem**.

All environment setup, data preparation, workflow design, execution, and interpretation were performed independently by me.  
Script implementation was supported by **automated/AI-assisted tools**, based on my instructions and scientific prompts.

---

##  Project Structure

```text
Virtual_Screening_Assignment/
│
├── Problem1_GoldenTriangle/
│   ├── data/
│   │   └── golden_triangle_infantshaliny.csv
│   ├── scripts/
│   │   └── problem1_filtering.py
│   └── results/
│       ├── golden_triangle_Filtered.csv
│       └── golden_triangle_infantshaliny.html
│
├── Problem2_ChemicalSpace/
│   ├── data/
│   │   └── chemical_space_infantshaliny.csv
│   ├── scripts/
│   │   └── problem2_plot.py
│   └── results/
│       ├── chemical_space_Descriptors.csv
│       └── chemical_space_infantshaliny_Plot.html
│
├── Problem3_RuleOfFive/
│   ├── data/
│   │   └── fda_drugs_infantshaliny.csv
│   ├── scripts/
│   │   └── problem3_lipinski.py
│   └── results/
│       ├── fda_drugs_infantshaliny_RO5_Analysis.csv
│       └── fda_drugs_infantshaliny_RO5_Violations_Interactive.html
│
└── report/
```
## Environment Setup
## To create and activate the conda environment:
```
conda env create -f chemenv_environment.yml
conda activate chemenv
conda install notebook -y
```
## How to Run the Workflow
## Golden Triangle Virtual Screening
```
cd Problem1_GoldenTriangle/scripts
python problem1_filtering.py
```
Input: ../data/golden_triangle_infantshaliny.csv

Output: Filtered CSV + HTML summary in ../results/
## Chemical Space Visualization
```
cd Problem2_ChemicalSpace/scripts
python problem2_plot.py
```
Generates a 2D scatter plot (logP vs TPSA)

Output: Interactive Plotly HTML file in ../results/

## Lipinski’s Rule of Five Analysis
```
cd Problem3_RuleOfFive/scripts
python problem3_lipinski.py
```
Input: ../data/fda_drugs_infantshaliny.csv

Output: RO5 violation summary (CSV + HTML) in ../results/
## Methods Summary
## Golden Triangle Filters

Criteria used:

Molecular Weight (MW) ≤ 350

logP ≤ 3.5

TPSA ≤ 90 Å²

0 RO5 violations

Compounds meeting these are considered lead-like.

## Chemical Space Mapping

Descriptors:

logP

TPSA

Visualization highlights:

Drug-like region

Outliers

Permeability & solubility balance

## Lipinski’s Rule of Five (RO5)

Evaluated thresholds:

MW ≤ 500

logP ≤ 5

HBD ≤ 5

HBA ≤ 10

Violations were identified and summarized.

## Tools & Libraries

Python 3.10

RDKit – molecular descriptors & filtering

pandas – data manipulation

plotly – interactive visualizations

Jupyter Notebook – optional exploration

Conda – environment management

## Author

Infant Shaliny Arokia Lourdu Elango
Bioinformatics Student · Cheminformatics & Drug Discovery Enthusiast

This repository demonstrates a compact computational drug discovery workflow using real molecular data and automated scripting tools.

## License

This project is released under the MIT License.
See the LICENSE file for more details.
