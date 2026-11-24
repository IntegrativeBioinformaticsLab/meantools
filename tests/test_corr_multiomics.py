"""Exercise the actual correlation-to-SQL caller in the full MEANtools env."""
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

import mutual_ranks


@unittest.skipUnless(importlib.util.find_spec('rdkit'), 'requires full MEANtools environment (rdkit)')
class CorrelationCallerTests(unittest.TestCase):
    def test_unimodal_and_bipartite_and_preserve_database(self):
        script = Path(__file__).resolve().parents[1] / 'corrMultiomics.py'
        help_result = subprocess.run([sys.executable, str(script), '--help'], capture_output=True, text=True)
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertNotIn('--network-type', help_result.stdout)
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary)
            matrix = pd.DataFrame([[1., 2., 3., 4., 5., 6.], [1., 3., 2., 6., 4., 5.],
                                   [6., 4., 5., 2., 3., 1.]], index=['A', 'B', 'C'],
                                  columns=[f'S{i}' for i in range(6)])
            matrix.to_csv(out / 'input.csv')
            # Different biological layers can have identical feature names
            # and even values. Their distinct input files preserve identity.
            matrix.to_csv(out / 'other_layer.csv')
            for network_type, single_rate in [('unimodal', False), ('bipartite', False),
                                               ('unimodal', True), ('bipartite', True)]:
                database = out / f'{network_type}_{single_rate}.sqlite'
                quantitation_matrix = out / ('input.csv' if network_type == 'unimodal' else 'other_layer.csv')
                command = [sys.executable, str(script), '-ft', str(out / 'input.csv'),
                           '-qm', str(quantitation_matrix), '-dn', str(database), '-tn', 'example',
                           '-mr', '-t', '1', '-c', '0', '-mdr', '25', '-d', '25']
                if single_rate:
                    # The single-rate branch is a programmatic API path: the
                    # existing CLI always supplies a nonempty multi-rate list.
                    code = ('import sys; sys.path.insert(0, ' + repr(str(script.parent)) + '); '
                            'import corrMultiomics as c; options=c.get_args(); '
                            'options.multi_decay_rates=[]; c.main(options)')
                    command = [sys.executable, '-c', code] + command[2:]
                environment = os.environ.copy()
                java = Path(sys.prefix) / 'lib/jvm/bin/java'
                if java.exists():
                    environment['PATH'] = str(java.parent) + os.pathsep + environment['PATH']
                    command.append('-cl')
                result = subprocess.run(command, cwd=out, env=environment, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                with sqlite3.connect(database) as connection:
                    correlations = pd.read_sql_query('SELECT * FROM example', connection)
                    edges = pd.read_sql_query('SELECT * FROM example_MR_edges', connection)
                    suffix = '' if single_rate else '_DR_25'
                    weights = pd.read_sql_query(f'SELECT * FROM example_MR_weights{suffix}', connection)
                    provenance = json.loads(connection.execute('SELECT provenance_json FROM example_MR_provenance').fetchone()[0])
                    if java.exists():
                        modules = pd.read_sql_query(f'SELECT * FROM example_clone{suffix}', connection)
                        self.assertFalse(modules.empty)
                        if network_type == 'bipartite':
                            members = set(' '.join(modules.Members).split())
                            self.assertTrue(all(x.startswith(('gene:', 'metabolite:')) for x in members))
                            self.assertIn('gene:A', members)
                            self.assertIn('metabolite:A', members)
                expected = mutual_ranks.get_MR_from_correlations(correlations, network_type=network_type)
                pd.testing.assert_frame_equal(edges, expected)
                self.assertEqual(len(edges), 3 if network_type == 'unimodal' else 9)
                np.testing.assert_allclose(weights.weight, np.exp(-(weights.MR - 1) / 25))
                self.assertEqual(provenance['network_type'], network_type)
                self.assertNotIn('network_type', provenance['settings']['options'])
                original = database.read_bytes()
                repeated = subprocess.run(command, cwd=out, env=environment, capture_output=True, text=True)
                self.assertNotEqual(repeated.returncode, 0)
                self.assertIn('requires a new output database', repeated.stderr)
                self.assertEqual(database.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
