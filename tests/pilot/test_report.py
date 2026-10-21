import unittest

from ibp.pilot.report import merge_loci, span


class MergeLociTest(unittest.TestCase):
    def test_nested_calls_are_one_locus(self):
        loci = merge_loci(["100-900", "200-300", "400-500"], 10_000)
        self.assertEqual(loci, [(100, 900, 3)])

    def test_separate_calls_stay_separate(self):
        loci = merge_loci(["100-200", "5000-6000"], 10_000)
        self.assertEqual([(a, b) for a, b, _ in loci], [(100, 200), (5000, 6000)])

    def test_pieces_across_origin_join(self):
        # two pieces at both ends of the chromosome = one event across the origin
        loci = merge_loci(["9500-10000", "0-300"], 10_000)
        self.assertEqual(len(loci), 1)
        a, b, n = loci[0]
        self.assertEqual((a, b, n), (9500, 300, 2))
        self.assertEqual(span(a, b, 10_000), 800)

    def test_wrapped_interval(self):
        loci = merge_loci(["9800-200", "100-150"], 10_000)
        self.assertEqual(loci, [(9800, 200, 2)])


if __name__ == "__main__":
    unittest.main()
