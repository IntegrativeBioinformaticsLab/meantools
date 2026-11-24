"""Numerical and streamed-pipeline parity tests for the HPC runner."""
import importlib.util
import hashlib
import io
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import mutual_ranks

SCRIPT = Path(__file__).resolve().parents[1] / 'corrMultiomics_HPC.py'
spec = importlib.util.spec_from_file_location('hpc', SCRIPT)
hpc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpc)


def options(**kwargs):
    defaults = dict(method='pearson', mad_filter=False, remove_zeros=False,
                    threads=2, block_size=3, correlation_cutoff=0.2, edge_weight_cutoff=0.01,
                    multi_decay_rates=[10, 25], decay_rate=5, clusterone=False, mutual_rank=True,
                    plot=False, annotation='', scratch_dir=None, mr_storage='disk', mr_memory_gb=None,
                    sqlite_table_name='test', if_exists='fail', network_type='unimodal')
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def run_ranks(data, nodes, storage, network_type='unimodal', **settings):
    with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as scratch, sqlite3.connect(':memory:') as output:
        opt = options(network_type=network_type, mr_storage=storage, **settings)
        if storage == 'disk':
            scratch.execute('CREATE TABLE edges (a INTEGER,b INTEGER,correlation REAL)')
            scratch.executemany('INSERT INTO edges VALUES (?,?,?)', data)
            hpc.mutual_rank(scratch, output, opt, np.asarray(nodes), tmp)
        else:
            chunks = [np.array(data[i:i+7], dtype=hpc.EDGE_DTYPE) for i in range(0, len(data), 7)]
            hpc.ram_mutual_rank(chunks, output, opt, np.asarray(nodes), tmp)
            assert chunks == []
        tables = [name for name, in output.execute("SELECT name FROM sqlite_master WHERE type='table'")]
        return {name: pd.read_sql_query('SELECT * FROM "' + name + '"', output) for name in tables}


class NumericalTests(unittest.TestCase):
    def test_submission_wrapper_passes_storage_once(self):
        # Execute the real wrapper with site commands stubbed: no job is
        # submitted and no installation/input/output paths are accessed.
        with tempfile.TemporaryDirectory() as tmp:
            captured = Path(tmp)/'arguments'
            shell = '''
module() { :; }
mamba() { :; }
srun() { printf '%s\\0' "$@" > "$CAPTURED_ARGUMENTS"; }
source "$WRAPPER_PATH"
'''
            for mode, budget in [('ram',''),('disk','1')]:
                env = os.environ.copy()
                env.update(MR_STORAGE=mode, MR_MEMORY_GB=budget, SLURM_JOB_ID='test',
                           SLURM_JOB_NODELIST='localhost', SLURM_CPUS_PER_TASK='1',
                           SLURM_TMPDIR=tmp, CAPTURED_ARGUMENTS=str(captured),
                           WRAPPER_PATH=str(SCRIPT.with_name('corrMultiomics.sh')))
                run = subprocess.run(['bash','-c',shell],env=env,capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                arguments = captured.read_bytes().decode().rstrip('\0').split('\0')
                self.assertEqual(arguments.count('--mr-storage'),1)
                self.assertEqual(arguments[arguments.index('--mr-storage')+1],mode)
                self.assertEqual('--mr-memory-gb' in arguments,bool(budget))
                self.assertIn('MR_STORAGE='+mode,run.stdout)
                self.assertEqual(list(Path(tmp).glob('corrMO.*')),[])

    def test_storage_selection_rejects_conflicting_repeated_options(self):
        required = ['corrMultiomics_HPC.py', '-ft', 'x.csv', '-qm', 'x.csv',
                    '-dn', 'out.sqlite', '-tn', 'test']
        for first, second in [('ram','disk'), ('disk','ram')]:
            stderr = io.StringIO()
            with patch.object(sys, 'argv', required + ['--mr-storage',first,'--mr-storage',second]), \
                    patch('sys.stderr', stderr), self.assertRaises(SystemExit) as error:
                hpc.get_args()
            self.assertEqual(error.exception.code, 2)
            self.assertIn('Conflicting --mr-storage', stderr.getvalue())
        for extra, expected in [([], 'disk'), (['--mr-storage','ram'], 'ram'),
                                (['--mr-storage=ram','--mr-storage','ram'], 'ram')]:
            with patch.object(sys, 'argv', required + extra):
                parsed = hpc.get_args()
            self.assertEqual(parsed.mr_storage, expected)
            self.assertNotIn('_explicit_mr_storage', vars(parsed))

    def test_ram_mode_opens_only_the_output_database(self):
        rng = np.random.RandomState(49)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.csv'
            pd.DataFrame(rng.randn(4,12), index=list('ABCD')).to_csv(path)
            database = str(Path(tmp)/'output.sqlite')
            opt = options(feature_table=path, quantitation_matrix=path,
                          sqlite_db_name=database, mr_storage='ram', scratch_dir=tmp)
            original_connect = sqlite3.connect
            opened, statements = [], []
            def tracked_connect(name, *args, **kwargs):
                opened.append(str(name))
                if Path(str(name)).name == 'ranks.sqlite':
                    self.fail('RAM mode opened a disk ranking database')
                connection = original_connect(name, *args, **kwargs)
                connection.set_trace_callback(statements.append)
                return connection
            with patch.object(hpc.sqlite3,'connect',side_effect=tracked_connect), \
                    patch.object(hpc,'mutual_rank',side_effect=AssertionError('RAM switched to disk')):
                hpc.main(opt)
            self.assertEqual(opened,[database])
            self.assertFalse(any('PRAGMA synchronous=OFF' in sql or 'PRAGMA journal_mode=OFF' in sql
                                 for sql in statements))

    def test_methods_match_scipy(self):
        rng = np.random.RandomState(22)
        raw_x, raw_y = rng.rand(7, 13) + .1, rng.rand(11, 13) + .1
        for method in ['pearson', 'pearsonlog', 'spearman']:
            a = np.round(raw_x, 1) if method == 'spearman' else raw_x
            b = np.round(raw_y, 1) if method == 'spearman' else raw_y
            opt = options(method=method, correlation_cutoff=0)
            x, _ = hpc.prepare(pd.DataFrame(a), opt)
            y, _ = hpc.prepare(pd.DataFrame(b), opt)
            for _, _, i, j, r, p in hpc.blocks(x, y, opt):
                for k in range(len(i)):
                    left, right = a[i[k]], b[j[k]]
                    if method == 'pearsonlog':
                        left, right = np.log10(left), np.log10(right)
                    expected = stats.spearmanr(left, right) if method == 'spearman' else stats.pearsonr(left, right)
                    np.testing.assert_allclose([r[k], p[k]], expected, atol=1e-12)

    def test_input_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'matrix.csv'
            for sep in [',', ';']:
                path.write_text(f'id{sep}s1{sep}s2\na{sep}1{sep}2\n')
                self.assertEqual(hpc.read_matrix(path).shape, (1, 2))
            path.write_bytes(b'id,s1\na,\x82')
            with self.assertRaisesRegex(ValueError, 'UTF-8'):
                hpc.read_matrix(path)

    def test_filters_precede_shared_sample_selection(self):
        matrix = pd.DataFrame([[1., 2., 4., 0.], [2., 4., 8., 9.]],
                              index=['zero_outside_shared', 'valid'], columns=['s1','s2','s3','extra'])
        _, names = hpc.prepare(matrix, options(remove_zeros=True), ['s1','s2','s3'])
        self.assertEqual(list(names), ['valid'])
        # Without -r, log-Pearson must not add a zero filter on unmatched samples.
        _, names = hpc.prepare(matrix, options(method='pearsonlog'), ['s1','s2','s3'])
        self.assertEqual(list(names), list(matrix.index))
        # Undefined logs are omitted just as NaN correlations fail the standard cutoff.
        matrix.loc['invalid'] = [-1., 2., 3., 1.]
        _, names = hpc.prepare(matrix, options(method='pearsonlog'), ['s1','s2','s3'])
        self.assertNotIn('invalid', names)

    def test_corrected_minimal_input_and_invariances(self):
        upper = [(0, 1, .9), (0, 2, .8), (1, 2, .7)]
        expected = mutual_ranks.get_MR_from_correlations(pd.DataFrame(
            [('A','B',.9), ('A','C',.8), ('B','C',.7)],
            columns=['metabolite','gene','correlation']), network_type='unimodal')
        for storage in ['disk', 'ram']:
            for seed in range(4):
                data = upper + [(b,a,r) for a,b,r in upper] + [(a,a,1.) for a in range(3)] + upper
                np.random.RandomState(seed).shuffle(data)
                for records in [upper, data]:
                    actual = run_ranks(records, ['A','B','C'], storage)['test_MR_edges']
                    pd.testing.assert_frame_equal(actual, expected)

    def test_ties_signed_scores_and_randomized_backend_parity(self):
        rng = np.random.RandomState(52)
        nodes = [f'N{i:02}' for i in range(12)]
        upper = [(a,b,float(rng.choice([-.9,-.2,.2,.8,1.])))
                 for a in range(12) for b in range(a+1,12)]
        data = upper + [(b,a,r) for a,b,r in upper] + [(a,a,1.) for a in range(12)]
        rng.shuffle(data)
        reference = pd.DataFrame([(nodes[a],nodes[b],r) for a,b,r in data],
                                 columns=['metabolite','gene','correlation'])
        expected = mutual_ranks.get_MR_from_correlations(reference, network_type='unimodal')
        for storage in ['disk','ram']:
            actual = run_ranks(data, nodes, storage)['test_MR_edges']
            pd.testing.assert_frame_equal(actual, expected)

    def test_duplicate_conflicts_and_tolerance(self):
        for storage in ['disk','ram']:
            for kind in ['unimodal','bipartite']:
                conflict = [(0,1,.9), (1,0,.8)] if kind == 'unimodal' else [(0,1,.9),(0,1,.8)]
                with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
                    run_ranks(conflict, ['A','B'], storage, kind)
                actual = run_ranks([(0,1,.9), (0,1,.9+1e-11)], ['A','B'], storage, kind)
                self.assertEqual(actual['test_MR_edges'].iloc[0].MR, 1.)
            with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
                run_ranks([(0,1,0.),(0,1,.9e-12),(0,1,1.8e-12)], ['A','B'], storage)

    def test_bipartite_overlapping_names_and_distinct_orientations(self):
        data = [(0,2,.9),(0,3,.8),(1,2,.7),(1,3,.6)]
        reference = pd.DataFrame([('A','A',.9),('A','B',.8),('B','A',.7),('B','B',.6)],
                                 columns=['metabolite','gene','correlation'])
        expected = mutual_ranks.get_MR_from_correlations(reference, network_type='bipartite')
        for storage in ['disk','ram']:
            actual = run_ranks(data + data, ['A','B','A','B'], storage, 'bipartite')['test_MR_edges']
            pd.testing.assert_frame_equal(actual, expected)

    def test_empty_diagonal_and_single_edge(self):
        for data in [[],[(0,0,1.)],[(0,1,-.4)]]:
            disk = run_ranks(data, ['A','B'], 'disk')
            ram = run_ranks(data, ['A','B'], 'ram')
            self.assertEqual(set(disk),set(ram))
            for table in disk:
                pd.testing.assert_frame_equal(disk[table],ram[table])
            self.assertEqual(len(disk['test_MR_edges']), 1 if data == [(0,1,-.4)] else 0)

    def test_pipeline_backends_and_rerun_protection(self):
        rng = np.random.RandomState(4)
        with tempfile.TemporaryDirectory() as tmp:
            x, y = Path(tmp)/'x.csv', Path(tmp)/'y.csv'
            # Shared labels across different files must remain separate nodes.
            pd.DataFrame(rng.randn(6,20),index=['F','A','B','C','D','E']).to_csv(x)
            pd.DataFrame(rng.randn(6,20),index=['E','B','C','A','F','D']).to_csv(y)
            for storage in ['disk','ram']:
                for unimodal in [False,True]:
                    opt = options(feature_table=x, quantitation_matrix=x if unimodal else y,
                        sqlite_db_name=str(Path(tmp)/f'{storage}_{unimodal}.sqlite'), mr_storage=storage)
                    hpc.main(opt)
                    with sqlite3.connect(opt.sqlite_db_name) as conn:
                        df = pd.read_sql_query('SELECT * FROM test',conn)
                        actual = pd.read_sql_query('SELECT * FROM test_MR_edges',conn)
                        expected = mutual_ranks.get_MR_from_correlations(df, network_type='unimodal' if unimodal else 'bipartite')
                        pd.testing.assert_frame_equal(actual,expected)
                        provenance = json.loads(conn.execute('SELECT provenance_json FROM test_MR_provenance').fetchone()[0])
                        self.assertEqual(provenance['ranking'],mutual_ranks.RANKING_CONVENTIONS)
                        self.assertEqual(provenance['network_type'],'unimodal' if unimodal else 'bipartite')
                        self.assertEqual(provenance['input_sha256'][str(x.resolve())],hashlib.sha256(x.read_bytes()).hexdigest())
                        for dr in [10,25]:
                            weights = pd.read_sql_query(f'SELECT * FROM test_MR_weights_DR_{dr}',conn)
                            np.testing.assert_allclose(weights.weight,np.exp(-(weights.MR-1)/dr))
                    original = Path(opt.sqlite_db_name).read_bytes()
                    opt.if_exists = 'replace'
                    with self.assertRaisesRegex(FileExistsError,'requires a new output database'):
                        hpc.main(opt)
                    self.assertEqual(Path(opt.sqlite_db_name).read_bytes(),original)

    def test_ram_budget_and_spearman_skips_mutual_ranks(self):
        with self.assertRaisesRegex(ValueError, 'exceeding'):
            hpc.check_ram_budget(1000000,.2)
        rng = np.random.RandomState(16)
        with tempfile.TemporaryDirectory() as tmp:
            x = Path(tmp)/'x.csv'
            pd.DataFrame(rng.randn(4,10),index=list('ABCD')).to_csv(x)
            opt = options(feature_table=x,quantitation_matrix=x,sqlite_db_name=str(Path(tmp)/'out.sqlite'),
                          mr_storage='ram',mr_memory_gb=1.,method='spearman',clusterone=True)
            hpc.main(opt)  # Like the standard script, this does not require Java.
            with sqlite3.connect(opt.sqlite_db_name) as conn:
                tables = {name for name, in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                self.assertIn('test',tables)
                self.assertFalse(any('_MR_' in name or '_clone' in name for name in tables))
            self.assertFalse(list(Path(tmp).glob('**/ranks.sqlite')))

    def test_cluster_failure_and_preserved_labels(self):
        with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as conn:
            out = Path(tmp)/'out.csv'
            opt = options(clusterone_jar='jar',java_heap='1g')
            with patch.object(hpc.subprocess,'run',return_value=subprocess.CompletedProcess([],1,stderr='bad network')):
                with self.assertRaisesRegex(RuntimeError,'bad network'):
                    hpc.cluster(opt,'weights',out,conn,'clusters')
            def fake_run(cmd, stdout, **kwargs):
                stdout.write('Cluster,Members\n1,metabolite:A gene:A\n')
                return subprocess.CompletedProcess(cmd,0,stderr='')
            with patch.object(hpc.subprocess,'run',side_effect=fake_run):
                hpc.cluster(opt,'weights',out,conn,'clusters')
            self.assertEqual(conn.execute('SELECT Members FROM clusters').fetchone()[0],'metabolite:A gene:A')


@unittest.skipUnless(importlib.util.find_spec('rdkit'), 'requires full MEANtools environment')
class StandardRunnerParityTests(unittest.TestCase):
    def test_actual_cli_tables_and_modules_match_standard(self):
        rng = np.random.RandomState(81)
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            samples = [f's{i}' for i in range(12)]
            # Intentionally unsorted row IDs exercise lexical tie/canonical ordering.
            pd.DataFrame(rng.rand(4,12)+.1,index=list('DACB'),columns=samples).to_csv(folder/'features.csv')
            pd.DataFrame(rng.rand(4,12)+.1,index=list('BDAC'),columns=samples).to_csv(folder/'genes.csv')
            env = os.environ.copy()
            java = Path(sys.prefix)/'lib/jvm/bin/java'
            if java.exists():
                env['PATH'] = str(java.parent)+os.pathsep+env['PATH']
            for method, unimodal, single in [('pearson',True,False),('pearson',False,True),
                                             ('pearsonlog',False,False),('spearman',False,False)]:
                common = ['-ft',str(folder/'features.csv'),'-qm',str(folder/('features.csv' if unimodal else 'genes.csv')),
                          '-tn','test','-mr','-t','1','-c','0.05','-mdr','10','25','-d','25','-m',method]
                if java.exists():
                    common += ['-cl']
                if single:
                    def command(script, db, extra):
                        code = ('import sys; sys.path.insert(0,'+repr(str(script.parent))+'); '
                                'import '+script.stem+' as c; options=c.get_args(); '
                                'options.multi_decay_rates=[]; c.main(options)')
                        return [sys.executable,'-c',code]+common+['-dn',str(db)]+extra
                else:
                    def command(script, db, extra):
                        return [sys.executable,str(script)]+common+['-dn',str(db)]+extra
                baseline = folder/f'standard_{method}_{unimodal}.sqlite'
                run = subprocess.run(command(ROOT/'corrMultiomics.py',baseline,[]),cwd=folder,env=env,capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                for storage in ['disk','ram']:
                    database = folder/f'{storage}_{method}_{unimodal}.sqlite'
                    run = subprocess.run(command(SCRIPT,database,['--mr-storage',storage,'--block-size','2']),
                                         cwd=folder,env=env,capture_output=True,text=True)
                    self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                    with sqlite3.connect(baseline) as reference, sqlite3.connect(database) as actual:
                        tables = [name for name, in reference.execute("SELECT name FROM sqlite_master WHERE type='table'")]
                        for table in tables:
                            left = pd.read_sql_query('SELECT * FROM '+table,reference)
                            right = pd.read_sql_query('SELECT * FROM '+table,actual)
                            if table.endswith('_provenance'):
                                p, q = json.loads(left.iloc[0].provenance_json),json.loads(right.iloc[0].provenance_json)
                                self.assertEqual(p['ranking'],q['ranking'])
                                self.assertEqual(p['network_type'],q['network_type'])
                            elif '_clone' in table:
                                self.assertEqual({frozenset(m.split()) for m in left.Members},
                                                 {frozenset(m.split()) for m in right.Members})
                            else:
                                order = ['metabolite','gene']
                                left,right = [frame.sort_values(order).reset_index(drop=True) for frame in [left,right]]
                                pd.testing.assert_frame_equal(left,right,atol=1e-12,rtol=1e-12)
                        if method == 'spearman':
                            actual_tables = {name for name, in actual.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                            self.assertFalse(any('_MR_' in name or '_clone' in name for name in actual_tables))


if __name__ == '__main__':
    unittest.main()
