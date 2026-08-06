from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from experimental_extensions.run_manifold import run


class ManifoldRegressionTest(unittest.TestCase):
    def test_pt3_gainsnzn_regression(self):
        root = Path(__file__).resolve().parents[1]
        config = root / "pt3_gainsnzn_regression.json"
        with tempfile.TemporaryDirectory() as temp_dir:
            summary = run(config, Path(temp_dir))
        self.assertEqual(summary["n_compositions"], 165)
        self.assertEqual(summary["n_pair_parameters"], 6)
        self.assertEqual(summary["design_rank"], 6)
        self.assertAlmostEqual(summary["training_R2"], 0.9992773421, places=8)
        self.assertAlmostEqual(
            summary["nonendmember_LOOCV_metrics"]["RMSE"],
            0.1422504672,
            places=8,
        )
        self.assertEqual(summary["transferability_claim"], "none")


if __name__ == "__main__":
    unittest.main()

