# Python NetMHCpan Wrapper

**A clean and efficient Python wrapper for NetMHCpan-4.2**

Seamless integration of NetMHCpan binding predictions into pandas-based bioinformatics workflows.

## Description

This package provides a user-friendly Python interface to **NetMHCpan-4.2** (standalone version).  
It handles peptide-MHC class I binding affinity predictions while solving common real-world problems such as inconsistent HLA nomenclature and missing HLA information.

The wrapper supports two main prediction regimes: **single allele** and **pan/supergroup** prediction, with built-in parallel processing for speed.

## Goals

- Make NetMHCpan easy to use inside Python/pandas scripts and notebooks
- Automatically fix inconsistent HLA allele names (e.g. `HLA-A*01:01`, `A0101`, `A01`)
- Enable pan-prediction when only broad allele families are known

## Supported Regimes

### 1. Single Mode (`single`)
- Requires an HLA column
- Automatically standardizes allele names using a  mapping file
- Predicts binding for the correctly matched allele

**Best for**: Standard epitope-HLA tables from IEDB, VDJdb, or custom datasets.

### 2. Pan / Supergroup Mode (`pan`)
- Works with or without a helper (supergroup) column
- If a valid family is given (e.g. `HLA-A01`, `HLA-B15`, `HLA-E01`), searches only inside that family
- If no family or unknown family → performs a smart **two-round search**:
  1. First round: finds the best supergroup using representative alleles
  2. Second round: finds the best allele inside the winning supergroup
- Gracefully handles completely HLA-free epitope lists

**Best for**: Large epitope sets with partial or no HLA information.


## Directory Structure
```text
mhc_netmhcpan_wrapper/
├── src/
│   ├── __init__.py
│   ├── config.py
│   ├── run_netmhcpan.py          # Main high-level function
│   ├── single_predictor.py
│   └── pan_predictor.py
├── datasets/
│   ├── allele_nomenclature_mapping.txt
│   └── mhc_supergroups.txt
├── tmp/                          # Auto-created for temporary files
├── testing.ipynb                 # Usage examples
├── Netmchpan_wrapper_development_notebook.ipynb
├── README.md
└── ...
```

***

## Quick Start Tutorial

### 1. Configuration
Edit `src/config.py` and set the number of cores:
```python
from pathlib import Path

TMPDIR = Path('../tmp')
PATH_TO_SUPERGROUPS = Path('../datasets/mhc_supergroups.txt')
PATH_TO_MAPPING = Path('../datasets/allele_nomenclature_mapping.txt')
N_CORES = 8
```

### 2. Basic Usage
```python
import pandas as pd
from pathlib import Path
from src.run_netmhcpan import run_netmhcpan

PATH_TO_NETMHCPAN = Path("/path/to/netMHCpan-4.2bstatic.Linux/netMHCpan-4.2")

# === Single mode ===
df = pd.read_csv("your_data.tsv", sep="\t")

result = run_netmhcpan(
    prediction_mode='single',
    path_to_netmhcpan=str(PATH_TO_NETMHCPAN),
    df=df,
    epitope_colname='antigen.epitope',
    hla_column='hla'
)

# === Pan mode ===
pan_df = pd.DataFrame({
    'epitope': ['LIDGIFLRY', 'VMADRTRHL', 'ANADLEVKI'],
    'family': [None, 'HLA-E01', 'HLA-C06']
})

result_pan = run_netmhcpan(
    prediction_mode='pan',
    path_to_netmhcpan=str(PATH_TO_NETMHCPAN),
    df=pan_df,
    epitope_colname='epitope',
    supergroup_column='family'   # can be None for pure pan search
)
```

***

### Using Classes Directly (Advanced)
```python
from src.single_predictor import SinglePredictor

single = SinglePredictor(
    path_to_netmhcpan=PATH_TO_NETMHCPAN,
    path_to_mapping=Path("datasets/allele_nomenclature_mapping.txt"),
    tmpdir=Path("tmp"),
    n_cores=8
)

df = single.match_hla(df, 'hla')
df = single.predict_affinity_dataframe(df, 'epitope', 'matched_hla')
```

***

## Important Notes
- You must have **NetMHCpan-4.2** installed and working from the command line (`netMHCpan -h`)
- The wrapper calls the binary via `subprocess`
- Predictions are run in **BA mode** (`-BA`)
- For very large datasets, adjust `N_CORES` according to your available CPU and memory

***

## License
This wrapper is available under GNU General Public License.

NetMHCpan is distributed by DTU Health Tech (https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) — please respect their licensing terms.