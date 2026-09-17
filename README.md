# LENS

LENS uses partial experimental interaction labels to calibrate microbiome network scores and rank unmeasured taxon pairs.
The model combines PLNNetwork, Poisson GLMNet, and SparCC scores with abundance features through centred ridge logistic regression.

This repository contains the code, prepared data, and results for the Butyrate assembly, Carlström, and Schäfer systems.

## Main files

- [Model](scripts/lens.py)
- [Prepared data](cleaned_data/)
- [Results](results/)
- [Manuscript and figures](paper/)
- [Data sources](external_data/PROVENANCE.md)

## Quick start

Run from the repository root using Python 3.12:

```sh
python3 -m pip install -r requirements.txt
python3 scripts/check_reproduction.py
```

This checks the saved predictions without overwriting results.
Prepared network scores are included, so network refitting is optional.

See the [reproduction guide](docs/reproduction.md) for analysis and figure commands.
