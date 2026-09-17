#!/usr/bin/env python3
"""Reconstruct the published 80% predictions without overwriting results."""
import numpy as np
import pandas as pd

from lens import ROOT, SYSTEMS, FEATURES, FIXED_FALLBACK_FEATURES, load_system, fit_score, CANDIDATES
from export_final_sparse_pr_predictions import run_system


def main():
    saved = pd.read_csv(ROOT / "results/final_sparse_pr_predictions_80pct.csv")
    for system in SYSTEMS:
        frame = load_system(system)
        assert frame.pair_id.is_unique, system
        assert set(frame.interaction_label.unique()) == {0, 1}, system
        assert len(FEATURES) == 14
        for label in (0, 1):
            training = frame[frame.interaction_label == label].head(10)
            np.testing.assert_allclose(
                fit_score(CANDIDATES[0], training, frame.head(5), 123),
                frame.head(5)[FIXED_FALLBACK_FEATURES].mean(axis=1),
            )
        actual = pd.DataFrame(run_system(system)).sort_values("pair_id")
        expected = saved[saved.analysis_set == system].sort_values("pair_id")
        assert list(actual.pair_id) == list(expected.pair_id), system
        for column in ["fold", "selected_setting", "interaction_label"]:
            assert np.array_equal(actual[column], expected[column]), (system, column)
        np.testing.assert_allclose(actual.supervised_score, expected.supervised_score, atol=1e-10, rtol=0)
        print(f"PASS {system}: {len(actual)} predictions, unchanged folds and penalty choices", flush=True)


if __name__ == "__main__":
    main()
