# Reproduction guide

**Supervised Calibration of Unsupervised Microbiome Network Scores with Partial Experimental Measurements**

LENS uses a small set of measured interactions to update a ranking derived from microbial abundance data.
This repository contains the current binary model and the three systems used in the manuscript.
Start with the prepared data and saved results; downloading raw sequencing data or refitting network estimators is optional.

## Quick start

Python requires NumPy, pandas, SciPy, and scikit-learn.
The analysis was run with Python 3.12.
Install dependencies in your own environment, then run these commands from the repository root:

```sh
python3 -m pip install -r requirements.txt
python3 scripts/check_reproduction.py
```

The check reconstructs all 2,617 predictions at the 80% budget and compares them with the saved results.
It does not overwrite data or launch network refits.
If using the existing local package directory, prefix Python commands with `PYTHONPATH=python_library`.

Recreate Figures 3–6 from the saved results:

```sh
Rscript scripts/plot_unsupervised_full_auprc.R
Rscript scripts/plot_final_sparse_pr_curves.R
Rscript scripts/plot_butyrate_absolute_pair_budgets.R
Rscript scripts/plot_carlstrom_signed_network.R
```

Plotting requires R with data.table, ggplot2, and igraph.
The budget plotting filename is historical: it produces all three panels.
The network plotting filename refers to colours showing experimental signs; LENS itself predicts binary interaction evidence.

## Data

| System | Evaluated pairs | Supported interactions | Tested neutral pairs |
| --- | ---: | ---: | ---: |
| Butyrate assembly | 104 | 76 | 28 |
| Carlström | 989 | 123 | 866 |
| Schäfer | 1,524 | 38 | 1,486 |

The counts refer to tested pairs with matching abundance features, after exclusions.
Source truth tables can contain additional pairs or ambiguous outcomes.
Untested pairs are never treated as neutral.
The interaction class combines positive and negative experimental effects.

See [data provenance](../external_data/PROVENANCE.md) for sources and preparation scripts.
Original data remain subject to their source terms.

## Model

[scripts/lens.py](../scripts/lens.py) contains feature construction, loading, fitting, fallback, and penalty selection.

Inputs are three score percentiles—PLNNetwork, Poisson GLMNet, and SparCC—and 11 abundance features.
The features describe prevalence, joint detection, presence similarity, abundance associations, log ratio variance, and abundance contrasts.

LENS fits logistic regression with balanced class weights and a ridge penalty centred on equal weights for the three scores.
Reference weights for abundance features are zero.
Median imputation and standardization use training data only.
Three inner stratified folds select the penalty from 0.01, 0.1, 1, and 10 by average precision.
If either training class has fewer than three observations, the model returns the mean of the three score percentiles.
Outputs are used for ranking; probability calibration has not been established.

Five stratified outer folds reserve pair labels within each system.
Test pairs stay fixed across budgets, and measured training samples are nested.
The absolute budget experiment repeats the five folds 20 times.
Unsupervised estimators use abundance data without interaction labels and remain fixed across budgets.
These experiments assess prediction within a system, not transfer to a new biological system.

## Reproduce the analysis

Export the 80% predictions:

```sh
python3 scripts/export_final_sparse_pr_predictions.py
```

Run the longer budget experiments:

```sh
python3 scripts/evaluate_butyrate_absolute_pair_budgets.py
python3 scripts/evaluate_carlstrom_absolute_pair_budgets.py
python3 scripts/evaluate_schafer_absolute_pair_budgets.py
```

The scripts overwrite their corresponding result files.
The base seed is 20260821; preserved seed offsets keep the saved partitions unchanged.

Network scores are already supplied in `analysis_data/`.
To refit them, use `build_pair_features.R` with a system identifier and `run_onenet_benchmarks.R`.
Package installation is described by `scripts/install_r_packages.R`; version records are in `analysis_data/package_versions.csv` and `results/onenet_package_versions.csv`.
A fully pinned portable environment has not yet been supplied.

Use `Rscript scripts/install_r_packages.R` for the individual estimators.
Set `INSTALL_ONENET=true` before that command to request the recorded OneNet revision as well.
This optional installer has not been tested in a fresh environment.

For the external abundance resamples in Figure 3:

```sh
Rscript scripts/bootstrap_abundance_refits.R
ONENET_EXTERNAL_WORKERS=2 Rscript scripts/bootstrap_onenet_refits.R
```

These refits are expensive.
OneNet uses six available estimators, 30 internal resamples, mean frequency aggregation, and mean stability 0.8 for density alignment.
ZiLN was excluded following failures involving constant columns.
Five external OneNet refits succeeded for Butyrate and Carlström; only one succeeded for Schäfer.
Four Schäfer refits failed during Magma.
Their failure records remain in `results/onenet_bootstrap_5/`.
The Schäfer OneNet bar has no estimated SD; its terminal cap is a graphical endpoint.

## Repository map

| Location | Purpose |
| --- | --- |
| `scripts/lens.py` | Current LENS implementation |
| `scripts/check_reproduction.py` | Prediction regression check |
| `scripts/` | Current evaluation, plotting, and source preparation |
| `cleaned_data/` | Abundance and experimental truth CSVs |
| `analysis_data/` | Network scores, pair features, and network display predictions |
| `results/` | Current predictions, metrics, partitions, and fit failure records |
| `paper/` | Manuscript, bibliography, and six PNG figures |
| `external_data/PROVENANCE.md` | Source descriptions |
| `interaction_ground_truth/` | Butyrate source tables and extraction scripts |

Build the manuscript with `cd paper && latexmk -pdf main.tex`.
This does not update Overleaf.

Obsolete experiments have been removed from the working tree.
They remain recoverable from Git commit `b4bb7efb33ecf3634ef2d266ffce7591a514482e`; see [cleanup notes](release_notes.md).
Local dependencies, caches, downloads, and the separate research idea log are excluded from Git.
