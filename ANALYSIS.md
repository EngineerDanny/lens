# Pair Feature and Training Fold Preparation

## Retained methods

Five methods completed the same fitting procedure for all five studies:

| Method | Software | Main settings | Exported evidence |
|---|---|---|---|
| PLNNetwork | PLNmodels 1.2.2 | 30 penalties; EBIC selection | Partial correlation, percentile rank, selected edge |
| Poisson GLMNet | glmnet 4.1.8 | Leave one taxon out; five fold internal tuning; `lambda.1se`; mutual selection | Symmetric coefficient, percentile rank, selected edge |
| SPIEC-EASI | SpiecEasi 1.99.0 | Graphical lasso; 30 penalties; minimum ratio 0.05; 20 StARS resamples | Partial correlation, percentile rank, selected edge, stability |
| SPRING | SPRING 1.0.4 | Data specific penalty sequence; 20 penalties; 20 StARS resamples | Symmetric coefficient, percentile rank, selected edge, stability |
| SparCC | SpiecEasi 1.99.0 | 20 outer iterations; 10 inner iterations | Correlation, percentile rank, selected edge at absolute correlation 0.3 |

Package versions are recorded in `analysis_data/package_versions.csv`. The project uses `r_library` so these packages do not alter the system R library.

## Omitted methods

gCoda installed successfully through NetCoMi 1.2.0. Its PairInteraX fit was slow and repeatedly reached the internal maximum iteration limit. It was removed from every study to preserve a common feature schema.

HARMONIES, COZINE, ZiLN, FlashWeave, compositional graphical lasso, and OneNet were omitted from this first analysis. Their software requires research scripts, another language, a large dependency chain, or methods already represented among the retained estimators.

## Input preparation

Native abundance values remain in `cleaned_data`.

PLNNetwork, Poisson GLMNet, SPIEC-EASI, and SparCC received integer input. Integer abundance tables were used directly. Relative abundance and proportion tables were rescaled within each sample to a total of 10,000 and rounded. The generated pair tables record this choice in `count_conversion`.

SPRING received the native abundance values. Absolute abundance and CFU tables were marked quantitative. PairInteraX was also treated as quantitative because it has no zeros. The butyrate proportions used SPRING's compositional procedure.

## Pair features

Each file named `analysis_data/<study>_pair_features.csv` contains every possible undirected pair among the retained taxa. Scores were computed without experimental labels.

The primary predictor list contains 22 columns:

- Five score percentile columns
- Five selected edge indicators
- SPIEC-EASI stability
- SPRING stability
- Seven prevalence and joint observation features
- Number of samples
- Number of taxa
- Overall zero frequency

Raw scores and additional abundance summaries remain in the fold files for later checks. `training_data/model_columns.csv` identifies the primary predictors and response.

## Training folds

Each fold folder contains `train.csv` and `test.csv`. All rows from a held out study remain together.

| Fold | Held out study | Training rows | Test rows | Test positives | Test neutrals |
|---|---|---:|---:|---:|---:|
| fold1 | OMM12 | 1,828 | 66 | 66 | 0 |
| fold2 | OMM12 keystone 2023 | 1,867 | 27 | 25 | 2 |
| fold3 | PairInteraX | 207 | 1,687 | 1,441 | 246 |
| fold4 | Butyrate assembly 2021 | 1,790 | 104 | 76 | 28 |
| fold5 | Host fitness 2018 | 1,884 | 10 | 9 | 1 |
| fold6 | Both OMM studies | 1,801 | 93 | 91 | 2 |

The primary folds exclude 136 ambiguous PairInteraX pairs. Those pairs remain in `training_data/all_tested_pairs.csv`. The 1,894 rows with binary labels are stored in `training_data/all_labeled_pairs.csv`.

## Leakage controls

- Experimental labels were unavailable during network fitting and feature construction.
- A complete study is assigned to either training or testing in each fold.
- Score percentiles were calculated separately within each abundance study.
- Predictor scaling has not been applied. Any model must estimate scaling values from its training file and apply them unchanged to its test file.
- Untested pairs are absent from binary training and testing.

## Known warnings

SPRING reported that its supplied penalty range may be too small for OMM12 keystone. It selected no edges there, although its stability values and zero score remain recorded.

SPIEC-EASI reported a similar penalty range warning for host fitness. Its selected graph contains one edge with a near zero fitted partial correlation.

One PairInteraX Poisson GLMNet target returned an empty model after a convergence warning. The procedure retained zeros for that target and completed the remaining fits.

These results remain in the feature files because the same fixed procedure was used across studies. Sensitivity analysis should examine wider penalty ranges before final model fitting.

## CoNet analysis of the Schäfer perturbation data

CoNet 1.1.1.beta was run through its command line interface. The ensemble used Pearson, Spearman, Bray Curtis, Kullback Leibler, and mutual information. Taxa present in fewer than 20 samples were removed before fitting, following the occurrence filter used in the published CoNet workflow. This retained 53 taxa and 552 tested pairs, of which 23 are positive.

The initial graph included the 500 upper and 500 lower scoring edges for each measure. This value avoids overlap because 53 taxa provide 1,378 possible pairs. Edge significance used 100 row shuffle permutations with renormalization followed by 100 bootstrap samples. Brown's method combined measure specific p values, and Benjamini Hochberg correction was applied at 0.05. Experimental labels were used only after the final graph was produced.

The exact analysis is in `scripts/run_conet_schafer.sh`. Pair scores are in `analysis_data/schafer_phyllosphere_2022_conet_tested_pair_scores.csv`, and the comparison restricted to CoNet's 552 eligible pairs is in `results/schafer_phyllosphere_2022_network_performance_conet_universe.csv`.

CoNet selected 437 tested pairs. It recovered 20 of the 23 positives, with precision 0.046, recall 0.870, and F1 0.087. Its average precision was 0.037, below the positive prevalence of 0.042. The result indicates a dense selected graph with weak discrimination on this perturbation benchmark. CoNet p values underflowed to zero for most retained edges at 100 resamples, so its full ranking uses measure support with the corrected p value as a tie breaker and should be treated cautiously.

## Compact calibration experiments

Three smaller models were trained using the same five training systems and evaluated on the untouched Schäfer labels. The score model used PLNNetwork, Poisson GLMNet, and SparCC score percentiles. The observability model used minimum prevalence and joint prevalence. The compact combined model used the three scores, method agreement, minimum prevalence, and joint prevalence.

The score model achieved AUPRC 0.254 and AUROC 0.899 on all 1,524 Schäfer pairs. It recovered 11 positives among the 38 highest ranked pairs, corresponding to density matched F1 0.289. A simple mean of the three score percentiles achieved AUPRC 0.229, while SparCC alone achieved 0.127. Supervised weighting therefore gave a modest ranking improvement beyond untrained score averaging.

The observability model achieved AUPRC 0.018. The compact combined model achieved AUPRC 0.085. Its weaker result suggests that prevalence and native edge agreement did not transfer reliably to this perturbation experiment.

The score model assigned high probabilities to nearly every Schäfer pair, and its training-derived threshold selected all 1,524 pairs. Its ranking is useful, but its probabilities are not calibrated for the sparse Schäfer truth. The current evidence supports a ranking claim and does not support automatic graph selection from transferred probabilities.

The experiment is implemented in `scripts/run_compact_calibration_experiments.py`. Results are in `results/compact_calibration_experiments_schafer.csv`, with coefficients, tuning records, metadata, and pair predictions stored beside the other analysis outputs.

The compact predictor sets were specified after the initial Schäfer result had been examined. Schäfer labels did not enter numerical fitting or tuning, but the experiment is still exploratory because model development was informed by performance on that system. The three-score model must be frozen and evaluated on another experimental system before making a confirmatory superiority claim.
