# LENS repository instructions

## Current study

All work belongs in this directory.
The current model is binary centred ridge LENS in `scripts/lens.py`.
The manuscript evaluates partial experimental measurements within Butyrate assembly, Carlström, and Schäfer.
Earlier transfer, taxon panel, signed classification, and model search experiments were removed from this working tree.
Do not restore them as the primary design.

## Scientific and implementation rules

- Preserve tested neutral, supported interaction, ambiguous, and untested outcomes distinctly.
- Never use interaction labels to fit or tune unsupervised estimators.
- Keep outer test labels out of preprocessing, penalty selection, and training pair selection.
- Preserve the fixed test pairs and nested training samples across budgets.
- Fit imputation and scaling on training data only.
- Preserve seeds, including the recorded system offsets.
- Treat folds as pair partitions within a biological system, not independent experiments.
- The network display uses binary predictions; colours annotate experimental signs.
- Keep fit failure records and source provenance.
- Do not change saved scientific results during a refactor.
- Run `scripts/check_reproduction.py` after model refactoring.
- Write scientific figure code in R and keep manuscript figure assets as PNG.
- Do not edit the separate PLN prediction study or synchronize Overleaf unless requested.

## Writing

Use direct academic prose.
Avoid unnecessary compound words, “rather than,” “not X, but Y,” “context,” and constructions ending in “-level.”
Keep one complete sentence per LaTeX source line.
Keep captions short and explain figures in the main text.
When editing manuscript prose, use the humanizer skill and preserve citations and LaTeX syntax.

## Native file-picker fallback

If a picker fails, open a fresh picker, navigate to the verified containing folder, focus the file list, type the filename, verify the selection, and press Return.
Confirm the web upload field displays the filename before claiming upload success.
Do not bypass file permissions.
