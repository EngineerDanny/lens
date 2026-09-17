# LENS

**Supervised Calibration of Unsupervised Microbiome Network Scores with Partial Experimental Measurements**

LENS uses partial experimental interaction measurements to update a ranking derived from microbiome abundance data. The current study asks how prediction of unmeasured pairs changes as more interactions are measured within the same microbial system.

The manuscript is in [paper/main.tex](paper/main.tex). The earlier README is preserved in [docs/analysis_history.md](docs/analysis_history.md). Older scripts and instruction documents contain superseded plans; the current analysis is described below.

## Experimental systems

The evaluation includes experimentally tested pairs with matching abundance features.

| System | Tested pairs | Supported interactions | Tested neutral pairs | Experimental evidence |
| --- | ---: | ---: | ---: | --- |
| Butyrate assembly | 104 | 76 | 28 | Effects on abundance during community assembly |
| Carlström | 989 | 123 | 866 | Abundance changes after strain removal |
| Schäfer | 1,524 | 38 | 1,486 | Abundance changes after strain addition |

Untested pairs are excluded from binary evaluation. Supported interactions combine positive and negative effects. Original signs and directions are retained where available. A neutral label means the experiment did not meet its criterion for an interaction under the conditions tested.

## Model and evaluation

LENS combines score percentiles from PLNNetwork, Poisson GLMNet, and SparCC with 11 abundance features describing prevalence, joint detection, presence similarity, abundance associations, log ratio variance, and abundance contrasts.

The model is logistic regression with balanced class weights and a ridge penalty centred on equal weights for the three network scores. Reference coefficients for abundance features are zero. Training data determine median imputation and standardization; reference coefficients are adjusted to preserve the ranking after scaling.

Inner validation selects the penalty from `0.01`, `0.1`, `1`, and `10` using average precision. When either training class contains fewer than three observations, prediction falls back to the equal mean of the three score percentiles. This check also applies inside inner validation. Outputs are ranking scores; calibrated interaction probabilities have not been established.

The main comparison uses five stratified folds of pair identifiers within each system. Each pair receives a prediction from a model trained without its label. Within a fold, test pairs stay fixed across budgets and training samples are nested. Unsupervised network scores use abundance data without interaction labels and remain fixed across budgets.

The absolute budget experiment uses 20 repetitions of five folds. Budgets begin at 10 measured pairs and increase to all eligible training pairs. The base seed is `20260821`, with derived seeds recorded in the scripts. Average precision, reported as AUPRC, is primary; AUROC is secondary.

These experiments assess prediction within the same system. They do not establish transfer to independent biological systems. Shared taxa and overlapping training sets mean that folds and repetitions are not independent biological replicates.

## Current figures and analyses

- **Figure 1:** OneNet and LENS workflows.
- **Figure 2:** Partial measurements and the pair validation design.
- **Figure 3:** Unsupervised AUPRC across external abundance resamples.
- **Figure 4:** Precision–recall curves at the 80% budget.
- **Figure 5:** Performance across absolute measurement budgets.
- **Figure 6:** Carlström interaction selections annotated with experimental signs.

OneNet uses PLNnetwork, SPIEC-EASI, gCoda, EMtree, Magma, and SPRING, with 30 internal resamples, mean frequency aggregation, and a target mean stability of 0.8 for density alignment. ZiLN was excluded after failures involving constant columns. This is a documented adaptation of OneNet.

For Figure 3, external resamples use seed `20260825` and the same sample indices for all methods. Constant taxa are omitted from fitting where necessary; their pairs receive zero scores to retain the original scoring universe.

| System | Successful OneNet refits | Mean AUPRC | SD |
| --- | ---: | ---: | ---: |
| Butyrate | 5 | 0.735 | 0.057 |
| Carlström | 5 | 0.122 | 0.004 |
| Schäfer | 1 | 0.113 | Unavailable |

Four Schäfer refits failed during Magma with constant taxon columns. Its OneNet bar therefore shows one successful refit, conditional on fitting success. The terminal cap is a graphical endpoint, not a zero SD estimate. Other capped error bars show one SD across successful refits. Failures remain recorded in the CSV results.

Figure 6 uses binary LENS predictions at the 80% budget. Both methods select 125 pairs across the five test folds, using selection counts derived from training interaction frequencies. LENS recovers 40 supported interactions and OneNet recovers 11. Colours show experimental effects, not predicted signs.

The separate [density diagnostic](results/carlstrom_density_recovery/README.md) permits different selection counts. It uses a training F1 threshold for LENS and a fixed mean frequency greater than 0.9 for OneNet. This exploratory diagnostic is not the selection procedure in Figure 6.

The taxon panel analysis uses an earlier fixed L2 model and a different budget denominator. It does not directly evaluate final LENS. Earlier signed classification results also must not be attributed to the current binary model.

## Repository layout

| Directory | Contents |
| --- | --- |
| `paper/` | LaTeX source, bibliography, and PNG figures |
| `scripts/` | Fitting, evaluation, and R plotting scripts |
| `cleaned_data/` | Prepared abundance tables and experimental labels |
| `analysis_data/` | Pair features, network scores, and predictions |
| `results/` | Metrics, audits, metadata, and diagnostic results |
| `docs/` | Historical analysis notes |
| `research_ideas/` | Literature search and research idea log |
| `external_data/`, `interaction_ground_truth/` | Local source material and provenance records |
| `r_library/`, `python_library/` | Local dependencies, excluded from Git |

## Reproduction

Run commands from the repository root. Python requires NumPy, pandas, SciPy, and scikit-learn. The existing local environment works with Python 3.12. R requires ggplot2, data.table, igraph, and the network packages used by each fitting script. Package records include `results/onenet_package_versions.csv`; a fully pinned portable environment has not yet been provided.

The main model code is in `scripts/optimize_supervised_model.py`. Penalty selection is in `scripts/evaluate_optimized_sparse_logistic.py`; its filename is historical, and its current candidates are centred ridge models. Abundance feature definitions are in `scripts/evaluate_direct_pair_features.py`.

With prepared data and network scores available:

```sh
PYTHONPATH=python_library python3.12 scripts/export_final_sparse_pr_predictions.py
PYTHONPATH=python_library python3.12 scripts/evaluate_butyrate_absolute_pair_budgets.py
PYTHONPATH=python_library python3.12 scripts/evaluate_carlstrom_absolute_pair_budgets.py
PYTHONPATH=python_library python3.12 scripts/evaluate_schafer_absolute_pair_budgets.py
```

Recreate result figures using the existing outputs:

```sh
Rscript scripts/plot_unsupervised_full_auprc.R
Rscript scripts/plot_final_sparse_pr_curves.R
Rscript scripts/plot_butyrate_absolute_pair_budgets.R
Rscript scripts/plot_carlstrom_signed_network.R
```

External OneNet refits are expensive. The runner reuses successful metrics and retries failed fits. Schäfer retries can reproduce the documented Magma error. Replotting existing results does not require refitting.

```sh
ONENET_EXTERNAL_WORKERS=2 Rscript scripts/bootstrap_onenet_refits.R
```

Build the local manuscript with a LaTeX installation providing `latexmk`:

```sh
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

These commands do not synchronize Overleaf. Rewrites discussed in chat may require separate application to the manuscript source.

## Git tracking

Track source scripts, manuscript text, bibliography, figure sources, publication PNGs, prepared CSV data, compact results, and provenance notes. The `.gitignore` excludes local dependencies, computational caches, binary fit objects, logs, temporary files, private environment files, and LaTeX build products.

Ignore rules do not remove files already committed. This update does not rewrite history or remove tracked data. Review `git status` and `git diff` before staging a research checkpoint. Third-party datasets remain subject to their original terms.
