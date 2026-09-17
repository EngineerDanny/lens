# Sources for the current study

Prepared CSV files are included in `cleaned_data/`.
Raw downloads are excluded from Git; analysis can begin with the prepared tables and saved network scores.
Source terms continue to apply.

## Butyrate assembly

- Article: https://doi.org/10.1038/s41467-021-22938-y
- Source release: `DesignSyntheticGutMicrobiomeAssemblyFunction-v1.0.zip`, from the Zenodo release identified by the article.
- Retained extracted sources: `interaction_ground_truth/butyrate_assembly_2021/processed/`.
- Preparation: `extract_butyrate_assembly_2021.R`, `build_truth_targets.R`, and `scripts/clean_interaction_data.py`.
- Abundance: experimentally measured proportions from assembled communities.
- Binary rule used here: an absolute aggregated log2 abundance ratio of at least 0.5.
- The cleaned tables contain 104 tested pairs, including 76 supported interactions.

## Carlström removal experiment

- Article: https://doi.org/10.1038/s41559-019-0994-z
- Source scripts and metadata: https://github.com/cmfield/carlstrom2019
- Sequencing deposit: ENA PRJEB32997.
- Preparation: `download_carlstrom_2019.py`, `process_carlstrom_2019_reads.py` or `process_carlstrom_2019_bulk.py`, `prepare_carlstrom_2019.py`, and `build_carlstrom_2019_truth.R`.
- Abundance scores use the intact community samples from experiment 53.
- Labels use the removal experiments, with DESeq2 adjusted p < 0.05 and an undirected interaction if either available direction meets that rule.
- The source truth table has 1,149 pairs. Matching the prepared abundance taxa retains 989 pairs, including 123 supported interactions.
- Effect signs describe observed abundance responses after removal; they are not predictions of facilitation or inhibition by LENS.

## Schäfer et al. Arabidopsis phyllosphere screen, 2022

- Study: *Mapping phyllosphere microbiota interactions in planta to establish genotype-phenotype relationships*
- Article: https://doi.org/10.1038/s41564-022-01132-w
- Downloaded: 2026-08-20
- Local source: `schafer_phyllosphere_2022/source`
- Abundance evidence: processed 16S amplicon counts and sample metadata from a 15-strain focal community challenged with individual added strains.
- Interaction evidence: adjusted Wald tests comparing each strain-addition condition with the unperturbed focal community. The published rule defines an effect at adjusted p <= 0.01.
- Use: addition experiment benchmark. Abundance profiles include the designed strain additions, so results must remain separate from observational abundance benchmarks.
- Harmonization: added strains are mapped to their reported ASV. Self-pairs and pairs with conflicting evidence after ASV mapping are excluded.

Key file checksums:

- `supplementary_data_1.xlsx`: `dc8cbb1d994f90cda54a434ae34575badb95394edf544a9a3c4e2166b18eefad`
- `supplementary_tables.xlsx`: `ead15eff3de5eef17fe2162671aee567a706c4bdd3a964dec3695341b5595da9`
- `source_data_figure_2.xlsx`: `cf832b2bcb8af85b9738289151568fdaed28efe0baf5a4c970f0c72b08e00867`
- `source_data_extended_figure_1.xlsx`: `98bf07b99ed689c4bc83cca3f3a8ebd3ae39a17dbf0506201b1d9f2c270f6efa`

## Reproduction boundaries

The source preparation scripts require their original downloads and, for Carlström, sequencing tools and DESeq2.
These are optional for reproducing the published analysis from the included prepared data.
The source manifest reports truth table counts before matching pair features; the main README reports the evaluated pairs.
Earlier source audits and excluded candidates remain in Git commit `b4bb7efb33ecf3634ef2d266ffce7591a514482e`.
