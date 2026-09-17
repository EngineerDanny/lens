# Carlström density recovery diagnostic

This exploratory experiment evaluates the number and identity of predicted interactions among the 989 experimentally tested Carlström pairs at the 80% measurement budget. It uses the existing five stratified outer folds and binary centred ridge LENS scores. Experimental truth contains 123 interactions (12.44%). Untested pairs are outside the evaluation.

Within each outer training set, three folds generate LENS predictions for training pairs whose labels are excluded from their fitted models. Each fit selects its own ridge penalty by inner AUPRC. The LENS threshold maximizes F1 on these training predictions; equal F1 values favour the higher threshold. The final threshold is applied to the existing outer test scores, selecting scores greater than or equal to the threshold.

OneNet selects mean selection frequencies strictly greater than 0.9, following the high frequency selection rule in Champion et al. (2024), https://doi.org/10.1371/journal.pcbi.1012627. Neither its scores nor this threshold use interaction labels. The stored OneNet run uses six estimators, 30 resamples, and an alignment stability target of 0.8. That alignment target is distinct from the final selection threshold. This is the published selection principle applied to our existing implementation, not an exact replication of every setting in the original paper.

The current rules were agreed before this revised diagnostic was run, following an earlier exploratory run that also tuned OneNet's threshold with training labels. The revised OneNet rule restores a fully unsupervised comparator. F1 evaluates interaction identification and does not directly optimize density error. This experiment measures the density induced by the stated selection rules; it does not establish the best possible density estimator.

LENS selected 250 pairs, including 55 supported interactions, giving precision 0.220, recall 0.447, and F1 0.295. OneNet selected one tested pair, which was neutral, giving precision, recall, and F1 of zero. Predicted densities were 25.28% and 0.10%, respectively, versus the observed 12.44%. Absolute pooled density errors were 12.84 and 12.34 percentage points. LENS overestimated density and OneNet underestimated it.

The mean absolute density error across folds was 13.67 percentage points for LENS and 12.34 for OneNet. Simply using the training interaction frequency to estimate test density gave an error of 0.28 percentage points. That baseline estimates a count without identifying which pairs interact. A small density error alone does not establish accurate interaction recovery. Precision is stored as missing for an individual fold when no pairs are selected; pooled OneNet precision is defined because one pair was selected.

The outer folds were stratified using labels, so their interaction frequencies are deliberately similar. The training frequency baseline benefits directly from that design. A stronger density experiment should use a prespecified set of splits without label stratification and assess both density error and interaction identification. The present experiment is a diagnostic using the existing manuscript folds. It does not validate transfer to new systems or complete network density when untested pairs lack truth.

Run from the repository root:

```sh
PYTHONPATH=python_library /opt/homebrew/bin/python3.12 scripts/evaluate_carlstrom_density_recovery.py
```

The script saves thresholds, fold metrics, nested training predictions, test selections, and pooled metrics in this folder. Existing manuscript figures are unchanged.
