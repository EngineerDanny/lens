# Supervised Calibration of Microbiome Network Scores

## Working paper identity

Working title:

> **Supervised Calibration of Unsupervised Microbiome Network Scores with Partial Experimental Measurements**

The title and abstract should center the methodological contribution expressed by the phrase “supervised calibration of unsupervised scores.” The main conceptual figure should compare the unsupervised OneNet workflow with the proposed supervised calibration workflow. Both begin with abundance data and multiple network estimators. OneNet aggregates estimator stability without experimental labels. The proposed workflow uses experimental labels from training systems to learn how estimator scores, prevalence, and joint detection predict interactions in a held-out experimental system.

## Project objective

Develop and evaluate a supervised calibration method that predicts experimentally supported microbial interactions by combining continuous scores from unsupervised network estimators with prevalence and joint-detection information.

The central research question is:

> Can experimental microbial interactions be predicted more accurately when unsupervised network scores are interpreted in light of taxon prevalence and joint observation?

The project must remain separate from the current SAGMB manuscript and the PLN versus GLMNet count-prediction study. All work for this project belongs in the current `network_inference` directory.

## Scientific framing

Unsupervised network estimators produce different statistical objects and use different graph-selection rules. Their continuous scores may contain useful interaction information even when their selected graphs have poor recall or incompatible densities.

The proposed method treats each estimator as a source of evidence. A small supervised model learns how those sources relate to experimentally measured interactions. Prevalence and zero-pattern features describe how much information is available for each taxon pair.

Do not claim that the zero fraction of an abundance table determines the density of its ecological network. The working hypothesis is narrower:

> Taxon prevalence and joint detection modify the reliability and interpretation of network-inference scores.

This is distinct from OneNet. OneNet forms an unsupervised stability consensus. The proposed model uses experimental labels to learn how estimator scores should be combined.

## Experimental benchmarks

Begin with these five datasets:

1. `omm12`
2. `omm12_keystone_2023`
3. `pairinterax`
4. `butyrate_assembly_2021`
5. `host_fitness_2018`

The source repository containing the existing processed data and scripts is:

`/Users/newuser/Projects/Personal/pln_eval`

Expected processed-data location:

`data/interaction_ground_truth/<dataset>/processed/`

Expected abundance files:

- OMM12 and OMM12 keystone: `community_absabundance_in_vitro.tsv.gz`
- PairInteraX, butyrate assembly, and host fitness: `abundance_matrix.tsv.gz`

Existing undirected truth file:

`truth_undirected.tsv.gz`

Copy or link source data only after documenting provenance. Do not modify the source repository's data or existing SAGMB files.

## Unit of analysis

The primary modeling table contains one row per experimentally tested taxon pair within an experimental dataset.

Required identifiers:

- Dataset
- Experimental study
- Taxon 1
- Taxon 2
- Direction, when available
- Experimental condition or medium, when available
- Truth type
- Tested-pair status

Primary response:

- `1`: experimentally supported interaction
- `0`: experimentally tested neutral pair
- Missing or excluded: pair was not experimentally tested

Never treat untested pairs as negative observations.

Preserve the original directional and signed outcomes. Create an undirected binary target only as a documented derived outcome.

## Truth types

Record the biological meaning of each label. At minimum, distinguish:

- Pairwise coculture effects
- Community dropout or add-back effects
- Community assembly effects
- Host fitness effects
- Broad undirected interaction summaries

Do not pool these outcomes without retaining truth type. Initially report performance separately. A pooled model must include study or truth-type structure.

## Unsupervised score generators

Initial methods:

- PLNNetwork
- GLMNet Poisson neighborhood selection
- SPIEC-EASI
- SPRING
- gCoda
- SparCC

Later candidates:

- FlashWeave
- HARMONIES
- COZINE
- ZiLN
- OneNet as an ensemble comparator

For every method, export continuous evidence whenever possible:

- Signed edge weight
- Absolute edge weight
- Within-dataset edge rank or percentile
- Bootstrap selection frequency
- Mean and standard deviation of bootstrap weights
- Native-selection indicator
- Penalty value or path position at which the edge enters

Experimental labels must remain unavailable during fitting, tuning, bootstrap estimation, and score construction for each unsupervised method.

## Prevalence and observability features

Compute the following for each taxon pair:

- Prevalence of taxon 1
- Prevalence of taxon 2
- Minimum and maximum prevalence
- Joint prevalence
- Number of jointly positive samples
- Number of samples containing only taxon 1
- Number of samples containing only taxon 2
- Number of samples containing neither taxon
- Jaccard similarity of presence patterns
- Mean abundance among positive observations
- Abundance variance or dispersion
- Sample size
- Number of taxa in the abundance table
- Overall zero frequency, used as context and not as a direct graph-density rule

Define presence using the raw or appropriately documented count-like table. Do not silently change the detection threshold.

## Optional biological features

Add these only when identifiers and provenance are reliable:

- Phylogenetic distance
- Shared taxonomy
- Genome similarity
- Metabolic pathway overlap
- Metabolic complementarity
- Oxygen tolerance
- Growth characteristics
- Substrate preferences

Evaluate biological features as a separate extension so their contribution can be distinguished from abundance-derived evidence.

## Primary supervised model

Use a small, interpretable model first:

- Regularized logistic regression for binary interaction prediction
- Standardize continuous predictors using training data only
- Account for correlated scores through regularization
- Include study or truth-type effects when datasets are pooled
- Consider a hierarchical logistic model after the basic pipeline is validated

The primary model should estimate an interaction probability for every tested pair. Preserve probabilities and do not reduce results immediately to a selected graph.

Nonlinear methods such as gradient boosting may be included as sensitivity analyses. They should not be the only supervised model.

## Baselines and ablations

Compare the supervised model against:

- Every individual unsupervised estimator
- Best individual estimator selected using training data only
- Mean standardized score
- Mean rank across estimators
- Majority vote among native graphs
- OneNet, if it can be run reliably
- Prevalence-only model
- Joint-detection-only model
- Organism-feature-only model, if biological features are available

Required ablations:

1. Unsupervised scores only
2. Prevalence and zero-pattern features only
3. Unsupervised scores plus prevalence features
4. Unsupervised scores plus prevalence and organism features

These comparisons test whether zero-pattern information improves the interpretation of network scores.

## Validation and leakage controls

Random edge splits are not acceptable as the primary validation design. Edges within a dataset share taxa, conditions, batches, and experimental procedures.

Use leave-one-experimental-system-out validation:

1. Hold out one complete experimental system.
2. Fit unsupervised methods on its abundance table without labels.
3. Train the calibration model using only the remaining systems.
4. Fix preprocessing, coefficients, hyperparameters, and thresholds.
5. Apply the frozen model to the held-out system.
6. Reveal the held-out labels only for evaluation.

OMM12 and OMM12 keystone may share organisms and experimental provenance. Include a stricter fold that holds both out together.

Where data permit, also evaluate:

- Leave-one-organism-out
- Leave-one-medium-out
- Leave-one-study-out
- Leave-one-truth-type-out

All score scaling must avoid leakage. Prefer within-abundance-dataset percentile ranks for method-specific scores. Any cross-dataset normalization parameters must be estimated using training systems only.

## Class imbalance and negative labels

Report the number of tested positives, tested neutrals, and untested pairs for every dataset.

If a study reports only positive interactions and leaves other pairs untested, do not construct ordinary binary negatives. Use one of the following:

- Restrict evaluation to explicitly tested pairs
- Use positive-unlabeled methods
- Exclude the dataset from supervised binary training while retaining it for positive-recovery analyses

Use class weights or appropriate regularization when tested positives and tested neutrals are imbalanced. Do not use accuracy as a primary metric.

## Evaluation

Primary metrics:

- AUPRC
- Precision
- Recall
- F1
- Number of selected edges
- Graph density

Probability-quality metrics:

- Brier score
- Calibration intercept and slope
- Reliability plots

Secondary metrics:

- AUROC, interpreted cautiously under class imbalance
- Precision at fixed edge budgets
- Recall at fixed edge budgets
- Density-matched F1
- Sign agreement
- Directional performance when labels and predictions permit it

Report full rankings before selecting a threshold.

Any probability or score threshold must be chosen using training data only. Acceptable strategies include nested cross-validation across training systems, a prespecified probability threshold, or a prespecified edge budget. Never maximize held-out F1.

## Data audit required before modeling

For every benchmark:

1. Identify all taxa in the abundance matrix.
2. Identify every experimentally tested pair.
3. Distinguish tested neutral pairs from untested pairs.
4. Preserve direction and sign.
5. Document how undirected truth was derived.
6. Calculate the number of possible, tested, positive, neutral, and untested pairs.
7. Calculate truth density among all possible pairs and among tested pairs.
8. Record truth type, experimental endpoint, medium, host, and community design.
9. Confirm that taxon identifiers match between abundance and truth files.
10. Record any exclusions or ambiguous mappings.

Do not begin supervised performance comparisons until this audit is complete.

## Initial deliverables

1. A machine-readable benchmark manifest describing all five datasets.
2. An audit table of tested, positive, neutral, and untested pairs.
3. A standardized edge-score schema shared by all estimators.
4. Exported continuous scores from the existing PLNNetwork, GLMNet, and SPIEC-EASI fits.
5. Implementations for SPRING, gCoda, and SparCC.
6. A pair-feature table containing estimator, prevalence, and joint-detection features.
7. A leakage-safe leave-one-system-out evaluation script.
8. Baseline, ablation, discrimination, and calibration results.
9. A concise methods document recording every preprocessing and tuning decision.

## Reproducibility requirements

- Use fixed and recorded random seeds.
- Record package versions and system dependencies.
- Cache computationally expensive unsupervised fits.
- Keep raw inputs, intermediate scores, fitted calibration models, and evaluation outputs separate.
- Use deterministic taxon-pair ordering.
- Store configuration in machine-readable files.
- Log exclusions and failed fits.
- Add tests for taxon alignment, pair ordering, symmetry, label status, and leakage boundaries.

## Scientific claims to avoid

Do not claim that:

- Statistical association establishes ecological causation.
- Matrix sparsity determines ecological network sparsity.
- Untested pairs are biological negatives.
- Pairwise coculture effects are equivalent to community or host effects.
- Random edge cross-validation measures transfer to new microbial systems.
- A supervised model is unsupervised because its input scores came from unsupervised estimators.

## Writing style

Use academic, direct, and natural language.

Avoid:

- “rather than”
- “not X, but Y” constructions
- Repeated negative constructions
- Unnecessary colons
- Unnecessary hyphenated expressions
- The word “context” unless no precise alternative works
- Constructions ending in “-level,” including “dataset-level,” “benchmark-level,” and “community-level”
- Compound words when a simpler phrase expresses the same meaning
- Excessive formulas or numerical detail in prose

When drafting manuscript text or reviewer responses:

- Name the exact section changed.
- Quote manuscript text briefly and exactly.
- Use ellipses for shortened quotations.
- Preserve LaTeX syntax and citations.
- Pass final wording through the humanizer skill.

## Immediate next action

Audit the five experimental truth files and their matching abundance matrices. Produce a table showing possible pairs, tested pairs, positive interactions, tested neutrals, untested pairs, direction availability, sign availability, and truth type. Do not fit the supervised calibration model until the tested-pair status is established.
