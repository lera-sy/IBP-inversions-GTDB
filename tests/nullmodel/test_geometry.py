import unittest

import numpy as np

from ibp.nullmodel.geometry import position_class, run, simulate_species


class GeometryNullTest(unittest.TestCase):
    def test_shorter_arc_and_floor(self):
        rng = np.random.default_rng(1)
        size, cls = simulate_species(2e6, 500, "random_breakpoints", rng, 6000)
        self.assertEqual(len(size), 500)
        self.assertTrue((size <= 1e6).all() and (size >= 6000).all())
        self.assertTrue(set(cls) <= {"local", "spanning"})

    def test_position_class(self):
        cls = position_class(np.array([100.0, 400.0]), np.array([300.0, 600.0]), 1000.0)
        self.assertEqual(list(cls), ["local", "spanning"])

    def test_random_breakpoints_slope_near_one(self):
        t = run([1, 10], n_species=40, n_events=20, replicates=20, seed=3)
        row = t[(t.null == "random_breakpoints") & (t.subset == "all") & (t.fit == "events")]
        self.assertAlmostEqual(float(row.slope_mean.iloc[0]), 1.0, delta=0.1)

    def test_small_fixed_sizes_give_flat_slope(self):
        t = run([1, 10], n_species=40, n_events=20, replicates=20, seed=3)
        row = t[(t.size_distribution == "lognormal_median_30kb") & (t.subset == "all") & (t.fit == "events")]
        self.assertLess(abs(float(row.slope_mean.iloc[0])), 0.15)


if __name__ == "__main__":
    unittest.main()
