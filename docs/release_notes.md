# Repository cleanup

The current release retains the binary centred ridge LENS analysis for Butyrate, Carlström, and Schäfer.
Obsolete transfer experiments, model searches, signed classifiers, density diagnostics, CoNet experiments, unused datasets, and unused figures were removed from the working tree.
Source preparation and provenance for the retained systems, current results, and OneNet failure records remain.

Earlier work is recoverable from Git commit `b4bb7efb33ecf3634ef2d266ffce7591a514482e`.
For example, `git show b4bb7ef:path/to/file` reads an old file without changing the working tree.
The cleanup does not rewrite Git history.
The research idea log remains on the local computer and is excluded from the release.
Ignored package libraries, raw downloads, and fit caches were left in place.

The required model functions were moved into `scripts/lens.py` before removing the old experiment scripts.
All 2,617 predictions at the 80% budget were reproduced, with maximum absolute differences below 1e-10.
Partitions, penalty selections, feature values, and pair order were unchanged.
The manuscript no longer describes the discarded taxon panel or signed classifier experiments.
