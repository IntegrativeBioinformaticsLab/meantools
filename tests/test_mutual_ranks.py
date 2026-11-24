import math
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

import mutual_ranks as mr


def correlations(rows):
    return pd.DataFrame(rows, columns=['metabolite', 'gene', 'correlation'])


def calculate(df, network_type='unimodal'):
    return mr.get_MR_from_correlations(df, network_type=network_type)


class MutualRankTests(unittest.TestCase):
    def setUp(self):
        self.upper = correlations([('A', 'B', .9), ('A', 'C', .8), ('B', 'C', .7)])
        self.expected = pd.DataFrame([
            ('A', 'B', 1.), ('A', 'C', math.sqrt(2)), ('B', 'C', 2.)
        ], columns=['metabolite', 'gene', 'MR'])

    def test_network_type_uses_file_identity_not_identifiers_or_contents(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first, other = root / 'first.csv', root / 'other.csv'
            self.upper.to_csv(first, index=False)
            self.upper.to_csv(other, index=False)
            self.assertEqual(mr.network_type_from_inputs(first, first), 'unimodal')
            self.assertEqual(mr.network_type_from_inputs(first, root / '.' / 'first.csv'), 'unimodal')
            self.assertEqual(mr.network_type_from_inputs(first, other), 'bipartite')
            symlink, hardlink = root / 'symlink.csv', root / 'hardlink.csv'
            symlink.symlink_to(first)
            os.link(first, hardlink)
            for alias in (symlink, hardlink):
                self.assertEqual(mr.network_type_from_inputs(first, alias), 'unimodal')

    def test_minimal_and_symmetric(self):
        reverse = self.upper.rename(columns={'metabolite': 'gene', 'gene': 'metabolite'})
        assert_frame_equal(calculate(self.upper), self.expected)
        assert_frame_equal(calculate(pd.concat([self.upper, reverse])), self.expected)

    def test_diagonal_duplicate_and_shuffle_invariance(self):
        reverse = self.upper.rename(columns={'metabolite': 'gene', 'gene': 'metabolite'})
        diagonal = correlations([(x, x, 1.) for x in 'ABC'])
        noisy = pd.concat([self.upper, reverse, diagonal, self.upper, diagonal])
        for seed in range(5):
            assert_frame_equal(calculate(noisy.sample(frac=1, random_state=seed)), self.expected)

    def test_exact_ties_use_partner_identifier_ordinal_ranks(self):
        tied = correlations([('A', 'C', .9), ('B', 'C', .9), ('A', 'B', .9)])
        # A: B=1,C=2; B: A=1,C=2; C: A=1,B=2.
        for seed in range(5):
            assert_frame_equal(calculate(tied.sample(frac=1, random_state=seed)), self.expected)

    def test_signed_negative_scores_are_not_ranked_by_absolute_value(self):
        negative = correlations([('A', 'B', .2), ('A', 'C', -.8), ('B', 'C', -.9)])
        assert_frame_equal(calculate(negative), self.expected)

    def test_conflicting_duplicates_including_reverse(self):
        for row in [('A', 'B', .8), ('B', 'A', .8)]:
            with self.assertRaisesRegex(ValueError, 'Conflicting duplicate.*1 pair'):
                calculate(pd.concat([self.upper, correlations([row])]))

    def test_tolerance_and_deterministic_representative(self):
        close = pd.concat([self.upper, correlations([('B', 'A', .9 + 1e-11)])])
        cleaned = mr.clean_correlations(close, network_type='unimodal')
        self.assertEqual(cleaned.iloc[0].correlation, .9 + 1e-11)
        assert_frame_equal(calculate(close), calculate(close.sample(frac=1, random_state=42)))
        # Validate the full range, not just each score against a middle value.
        chain = correlations([('A', 'B', x) for x in [0., .9e-12, 1.8e-12]])
        with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
            calculate(chain)

    def test_bipartite_disjoint_layers(self):
        df = correlations([('M1', 'G1', .9), ('M1', 'G2', .8),
                           ('M2', 'G1', .7), ('M2', 'G2', .6)])
        expected = pd.DataFrame([('M1', 'G1', 1.), ('M1', 'G2', math.sqrt(2)),
                                 ('M2', 'G1', math.sqrt(2)), ('M2', 'G2', 2.)],
                                columns=['metabolite', 'gene', 'MR'])
        assert_frame_equal(calculate(df, 'bipartite'), expected)
        assert_frame_equal(calculate(pd.concat([df, df]).sample(frac=1), 'bipartite'), expected)

    def test_bipartite_overlapping_identifiers_and_opposite_layer_edges(self):
        df = correlations([('A', 'A', .9), ('A', 'B', .8),
                           ('B', 'A', .7), ('B', 'B', .6)])
        edges = calculate(df, 'bipartite')
        np.testing.assert_allclose(edges.MR, [1., math.sqrt(2), math.sqrt(2), 2.])
        self.assertEqual(len(edges), 4)
        weights = mr.get_weight_from_MR(edges, .01, 25)
        exported = mr.clusterone_edges(weights, network_type='bipartite')
        self.assertEqual(exported.iloc[0].metabolite, 'metabolite:A')
        self.assertEqual(exported.iloc[0].gene, 'gene:A')
        self.assertTrue(set(exported.metabolite).isdisjoint(exported.gene))

    def test_bipartite_conflicting_repeated_edge(self):
        with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
            calculate(correlations([('A', 'A', .9), ('A', 'A', .8)]), 'bipartite')

    def test_empty_or_only_diagonal(self):
        for df in [correlations([]), correlations([('A', 'A', 1.)])]:
            output = calculate(df)
            self.assertEqual(list(output.columns), ['metabolite', 'gene', 'MR'])
            self.assertTrue(output.empty)

    def test_invalid_inputs_and_explicit_type(self):
        with self.assertRaises(TypeError):
            mr.get_MR_from_correlations(self.upper)
        with self.assertRaisesRegex(ValueError, 'network_type'):
            calculate(self.upper, 'auto')
        for row in [('A', 'B', np.nan), ('A', 'B', np.inf), (None, 'B', .9)]:
            with self.assertRaises(ValueError):
                calculate(correlations([row]))
        with self.assertRaisesRegex(ValueError, 'missing input endpoints'):
            mr.get_MR_from_correlations(self.upper, ['A'], network_type='unimodal')

    def test_inputs_unmodified_and_legacy_node_list(self):
        original = self.upper.copy(deep=True)
        assert_frame_equal(mr.get_MR_from_correlations(self.upper, ['C', 'A', 'B', 'A'],
                                                     network_type='unimodal'), self.expected)
        assert_frame_equal(self.upper, original)

    def test_weight_formula_and_strict_cutoff_unchanged(self):
        edges = calculate(self.upper)
        expected = np.exp(-(edges.MR - 1) / 25)
        output = mr.get_weight_from_MR(edges, expected.iloc[1], 25)
        self.assertEqual(len(output), 1)
        self.assertEqual(output.iloc[0].weight, 1.)


if __name__ == '__main__':
    unittest.main()
