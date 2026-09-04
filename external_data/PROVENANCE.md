# Public experimental candidates

## Retained source

### Friedman, Higgins and Gore soil microcosm, 2017

- Study: *Community structure follows simple assembly rules in microbial microcosms*
- Article: https://doi.org/10.1038/s41559-017-0109
- Supporting network analysis: https://doi.org/10.1038/s41467-017-02090-2
- Data repository: https://doi.org/10.5281/zenodo.8176044
- Downloaded: 2026-08-21
- Local source: `friedman_microcosm_2017/source`
- Archive SHA256: `d12276048892c048a507c63dcbba98e368438a289ef94c8757777757009aec27`
- Abundance evidence: endpoint absolute abundances from all trio cultures and all seven- and eight-species cultures. Absolute abundance was calculated by the study authors from total optical density and species fractions measured by plating.
- Interaction evidence: endpoint abundance of each focal species in pair culture compared with its monoculture endpoint abundance.
- Separation rule: network estimators use only the higher-order cultures. Monoculture and pair-culture measurements are reserved for experimental labels.
- Label rule fixed before network fitting: Welch test for each of 56 directed effects followed by Benjamini-Hochberg correction across all effects. Supported effects have `q <= 0.05`; the remaining tested effects are neutral.
- Exclusions: three incomplete higher-order endpoints and four zero-total cultures.
- License: the archived repository contains an MIT license.

### Schäfer et al. Arabidopsis phyllosphere screen, 2022

- Study: *Mapping phyllosphere microbiota interactions in planta to establish genotype-phenotype relationships*
- Article: https://doi.org/10.1038/s41564-022-01132-w
- Downloaded: 2026-08-20
- Local source: `schafer_phyllosphere_2022/source`
- Abundance evidence: processed 16S amplicon counts and sample metadata from a 15-strain focal community challenged with individual added strains.
- Interaction evidence: adjusted Wald tests comparing each strain-addition condition with the unperturbed focal community. The published rule defines an effect at adjusted p <= 0.01.
- Planned use: perturbation benchmark. Abundance profiles include the designed strain additions, so results must remain separate from observational abundance benchmarks.
- Harmonization: added strains are mapped to their reported ASV. Self-pairs and pairs with conflicting evidence after ASV mapping are excluded.

Key file checksums:

- `supplementary_data_1.xlsx`: `dc8cbb1d994f90cda54a434ae34575badb95394edf544a9a3c4e2166b18eefad`
- `supplementary_tables.xlsx`: `ead15eff3de5eef17fe2162671aee567a706c4bdd3a964dec3695341b5595da9`
- `source_data_figure_2.xlsx`: `cf832b2bcb8af85b9738289151568fdaed28efe0baf5a4c970f0c72b08e00867`
- `source_data_extended_figure_1.xlsx`: `98bf07b99ed689c4bc83cca3f3a8ebd3ae39a17dbf0506201b1d9f2c270f6efa`

### Wortel laboratory gut synthetic community, 2026

- Study: *Environmentally mediated interactions predict community assembly and invasion success in a gut microbiota synthetic community*
- Article: https://doi.org/10.1128/msystems.00113-26
- Data repository: https://github.com/WortelLab/SynCom
- Repository commit: `8d0138613eed585b3f98641953a03c407eb3142b`
- Commit date: 2026-03-25
- Downloaded: 2026-08-20
- Local source: `wortel_syncom_2026/source/SynCom`
- License for the article: CC BY 4.0. The GitHub repository does not contain a separate license file, so redistribution terms for the repository files should be confirmed before publication.
- Abundance evidence: qPCR measurements from 45 two species cultures, 11 three or four species cultures, and separate invasion experiments.
- Interaction evidence: growth rate and final population density effects measured in conditioned medium. These measurements are separate from the qPCR coculture measurements.
- Planned use: derive one abundance table from the coculture biological replicates and derive tested pair labels from the conditioned medium experiments. Keep invasion measurements as a separate sensitivity analysis.
- Important restriction: do not treat every nonzero reported coefficient as a binary interaction. The binary rule must use replicate uncertainty and must be fixed before evaluating network scores.

Key file checksums:

- `QPCR data/Files/In_vitro_data_with_raw.csv`: `d3a48fb6b333ff1174a8c2ee487560f75c40e9143a4e5a4617b122a9fa500ffc`
- `Parameters/output/a_b_parameters.csv`: `dcd6320b5b2cc42264406c1b19788b14383d1a53e58d6178ca21f41b84bc37c7`
- `Parameters/output/epsilons_RI.csv`: `1b7e368fe77d4f5bc88b373efadcff5d2873d0337ca1b9c66c8a75dd4eecfb33`

## CoNet software

- Software: CoNet 1.1.1.beta
- Official release page: https://apps.cytoscape.org/apps/conet
- Executable JAR: https://apps.cytoscape.org/download/conet/1.1.1.beta
- Downloaded: 2026-08-20
- Local file: `conet/conet-1.1.1.beta.jar`
- SHA256: `28e03f1cfbe33ce81c43c2f38902e7f12ba210b814d2ed83c694825535b18846`
- Runtime used: Java 17.0.1

## Removed candidates


The local copies and temporary audit outputs for the following candidates were deleted on 2026-08-20 after audit. The exclusion reasons are retained here.

### Ratzke et al. 2020

- Study: *Strength of species interactions determines biodiversity and stability in microbial communities*
- Repository: https://doi.org/10.5061/dryad.vdncjsxq9
- Reason for removal: the downloaded processed files contain total optical density and pH, without taxon resolved community abundance. The separate raw sequencing deposit is about 9.24 GB and would require a new processing and taxon alignment pipeline.

### Chang et al. 2023

- Study: *Emergent coexistence in multispecies microbial communities*
- Repository: https://doi.org/10.5061/dryad.bnzs7h4gb
- Reason for removal: strain to ESV mapping collapses experimental strain pairs into repeated or self pairs, some definitive outcomes conflict after mapping, and coexistence is not an experimentally tested neutral interaction. It is unsuitable for the primary positive versus neutral target.
