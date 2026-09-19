# Placement of interaction measurements

This experiment compares the current binary LENS model under restricted and distributed measurement designs.
It is a separate analysis; existing manuscript results are unchanged.

For each system, 20 random partitions reserve approximately 40% of taxa that appear in the evaluated pair table.
Test pairs contain two reserved taxa.
Restricted training uses pairs containing no reserved taxon.
Distributed training can use any other tested pair, including pairs connecting reserved taxa to other taxa.
The exact test pairs are excluded from both training pools.
All methods retain the complete abundance table.

Both designs receive the same number of labels and use the same test pairs within a partition.
A common random priority orders eligible pairs, with each design using its eligible subsequence.
Training samples are nested across budgets.
Budgets are capped at the smallest restricted pool across the 20 partitions, using pair availability without outcomes.
Thus the included partitions and OneNet reference do not change with the budget.

Partitions are drawn once without consulting outcomes.
Empty tests and tests containing only one class are recorded and excluded from both metric summaries, without redrawing.
Results are conditional on the remaining evaluable partitions.
The audit reports those exclusions explicitly.

The existing LENS architecture, features, inner pair validation and fallback are unchanged.
Fewer than three training examples in either class triggers the mean of the three network scores.
Inner validation remains the current pair validation procedure; it has not been redesigned to optimize transfer to reserved taxa.
Training class counts and coverage of test taxa are recorded, since equal label budgets need not yield equal class counts or organism coverage.

The figure shows grouped boxplots of AUPRC across evaluable partitions (17 for Butyrate and 20 for each other system).
The centre line is the median, the box spans the first and third quartiles, and whiskers extend to observed values within 1.5 times the interquartile range.
Values beyond the whiskers appear as points.
Budgets are displayed as ordered categories with equal spacing, and panel scales differ.
These partitions reuse one biological system and are not independent biological replicates.
The boxes describe variation among partitions, not confidence intervals or variation across biological systems.
Paired differences are computed on identical test sets; box overlap is not used to decide whether a difference exists.
OneNet and individual fixed estimators are scored separately, without label tuning or oracle selection.

Run from the repository root:

```sh
PYTHONPATH=python_library python3.12 scripts/evaluate_measurement_placement.py
Rscript scripts/plot_measurement_placement.R
```

The plan files record exact taxon partitions, test pairs, and training orders.
Prediction files support independent metric checks.
The metadata records seeds and input hashes.
No new unsupervised fits are required.

## Initial results

At the largest budget shared by all partitions, mean AUPRC was:

| System | Measured pairs | Restricted LENS | Distributed LENS | OneNet | Fixed mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| Butyrate | 13 | 0.721 | 0.745 | 0.755 | 0.736 |
| Carlström | 275 | 0.211 | 0.232 | 0.125 | 0.137 |
| Schäfer | 357 | 0.271 | 0.323 | 0.213 | 0.296 |

Butyrate retained 17 of 20 partitions; two had only one test class and one had no test pairs.
The other systems retained all 20 partitions.
The small common budget and sparse test coverage limit the Butyrate comparison.

In Carlström, both designs exceeded OneNet in all 20 partitions at the largest budget.
Restricted LENS also exceeded the fixed mean in 18 of 20 partitions, supporting useful updating for taxa without interaction labels in this system.
Distributed LENS had a higher mean than restricted LENS but won in only 10 of 20 partitions.
Its average advantage therefore does not establish consistent superiority of distributed measurements.

In Schäfer, distributed LENS exceeded restricted LENS in 14 of 20 partitions at the largest budget.
Restricted LENS had a lower mean AUPRC than the fixed reference, despite exceeding OneNet on average.
The results suggest that experimental coverage can help, but do not establish a universal allocation rule.
At five measured pairs both designs use the fallback in every partition; these values do not demonstrate supervised learning.

All 2,648 method evaluations were checked using independent AP and AUROC calculations from saved predictions.
The checker also verifies separation of training and test pairs, the restriction on training taxa, and fixed test membership and reference metrics across budgets.
Run `scripts/check_measurement_placement.py` to repeat these checks.
