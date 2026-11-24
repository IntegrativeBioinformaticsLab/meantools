#!/usr/bin/env python
"""Bounded-memory, single-node SLURM correlation and mutual-rank pipeline."""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
import sqlite3
import tempfile
import subprocess
import zipfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
from urllib.parse import quote as escape_identifier

# Reuse the standard runner's ranking policy and automatic layer selection.
# Resolve imports relative to this file so batch jobs may run from any directory.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import mutual_ranks

import numpy as np
import pandas as pd
from scipy import stats, special
try:
    from threadpoolctl import threadpool_limits
except ImportError:
    raise SystemExit("Missing threadpoolctl. Install in the active environment: python -m pip install threadpoolctl")

class StorageModeAction(argparse.Action):
    """Do not silently override a requested backend with a wrapper's default."""
    def __call__(self, parser, namespace, value, option_string=None):
        previous = getattr(namespace, '_explicit_mr_storage', None)
        if previous is not None and previous != value:
            parser.error(f'Conflicting --mr-storage values: {previous} and {value}. '
                         'Supply the option once; check MR_STORAGE and the submission script arguments.')
        setattr(namespace, '_explicit_mr_storage', value)
        setattr(namespace, self.dest, value)


def get_args():
    parser = argparse.ArgumentParser(prog='corrMultiomics',
                                     description='Correlate multi-omics and convert correlation coeffecient into edge weights!',
                                     epilog='Contact - kumarsaurabh.singh@maastrichtuniversity.nl')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-ft','--feature_table', default='', action='store', required=True, help='Metabolomics-based feature table.')
    required.add_argument('-qm','--quantitation_matrix', default='', action='store', required=False, help='Normalized expression table from the RNAseq data.')
    required.add_argument('-dn', '--sqlite_db_name', default='', action='store', required=True, help='Provide a name for the database!')
    required.add_argument('-tn', '--sqlite_table_name', default='', action='store', required=True, help='Provide a name for the database table.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-m', '--method', default='pearson', type=str, required=False, choices=['spearman', 'pearson', 'pearsonlog'], help='Default is the pearson correlation method; pearsonlog uses log-transformed arrays, where arrays with atleast a zero are always removed first regardless of other options.')
    optional.add_argument('-mr', '--mutual_rank', default=False, required=False, action='store_true', help='MR has more module predictive power than pearson or pearson log.')
    optional.add_argument('-cl', '--clusterone', default=False, required=False, action='store_true', help='Clusterone creates modules based on MR and edge weights.')
    optional.add_argument('-mad', '--mad_filter', default=False, required=False, action='store_true', help='Removes arrays with a MAD of 0.')
    optional.add_argument('-r', '--remove_zeros', default=False, required=False, action='store_true', help='Removes arrays with at least one 0.')
    optional.add_argument('-c', '--correlation_cutoff', default=0.1, required=False, type=float, help='Minimum correlation coeffecient cut-off. Default: 0.3')
    optional.add_argument('-w', '--edge_weight_cutoff', default=0.01, required=False, type=float, help='Minimum MR-derived edge weight. Default: 0.01')
    optional.add_argument('-d', '--decay_rate', default=5, required=False, type=int, help='Decay rate of default: 25 would be used.')
    optional.add_argument('-mdr', '--multi_decay_rates', default=[5,10,25], required=False, nargs='+', type=int, help='Decay rate could be any value but generally it is either {5, 10, 25, 50, 100}. In this case all 5 values will be used separately to estimate edges. Usage:: -mdr 5 10 25 50 100')
    optional.add_argument('-a', '--annotation', default='', type=str, required=False, help='Comma-delimited two-columns file with annotations. No header.')
    optional.add_argument('-p', '--plot', default=False, required=False, action='store_true', help='Plots showing the correlation will be created.')
    optional.add_argument('-t', '--threads', default=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')), type=int, required=False)
    optional.add_argument('-v','--verbose', default=False, action='store_true', required=False)
    optional.add_argument('--block-size', type=int, default=128, help='Feature rows per correlation block.')
    optional.add_argument('--mr-storage', choices=['ram', 'disk'], default='disk', action=StorageModeAction,
                          help='Mutual-rank intermediate storage: NumPy RAM arrays or scratch SQLite (default).')
    optional.add_argument('--mr-memory-gb', type=float, default=None,
                          help='RAM-mode estimated working-memory budget in GiB; defaults to 60%% of SLURM_MEM_PER_NODE when set.')
    optional.add_argument('--scratch-dir', default=os.environ.get('SLURM_TMPDIR'), help='Local scratch for ranking and ClusterOne files.')
    optional.add_argument('--clusterone-jar', default=str(next((base / 'cluster_one-1.0.jar' for base in (Path(__file__).resolve().parent, Path(__file__).resolve().parents[1]) if (base / 'cluster_one-1.0.jar').is_file()), Path(__file__).resolve().parent / 'cluster_one-1.0.jar')))
    optional.add_argument('--java-heap', default='8g', help='ClusterOne Java maximum heap; leave room in the SLURM memory allocation.')
    optional.add_argument('--if-exists', choices=['fail', 'append', 'replace'], default='fail', help='Correlation-only table collisions; MR runs always require a new database.')
    options = parser.parse_args()
    vars(options).pop('_explicit_mr_storage', None)
    return options




def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + message, flush=True)


@contextmanager
def stage(name):
    start = time.monotonic()
    log(name + ' started')
    try:
        yield
    finally:
        log(f'{name} elapsed={time.monotonic() - start:.2f}s')


def quote(name):
    return '"' + name.replace('"', '""') + '"'


def file_sha256(path):
    """Hash large input matrices without loading a second full file into RAM."""
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_matrix(path):
    if zipfile.is_zipfile(path):
        raise ValueError(f'{path}: spreadsheet/archive rather than text CSV. Export as UTF-8 CSV.')
    try:
        df = pd.read_csv(path, index_col=0, sep=None, engine='python', encoding='utf-8-sig')
    except UnicodeDecodeError as exc:
        raise ValueError(f'{path}: export this file as UTF-8 CSV.') from exc
    except (pd.errors.ParserError, pd.errors.EmptyDataError, csv.Error) as exc:
        raise ValueError(f'{path}: empty or malformed CSV; check delimiter and quoting: {exc}') from exc
    if df.empty or not df.index.is_unique or df.index.hasnans:
        raise ValueError(f'{path}: require nonempty data and unique, nonmissing row IDs.')
    df.index = df.index.astype(str)
    if not df.index.is_unique:
        raise ValueError(f'{path}: IDs are not unique after converting them to strings.')
    return df


def prepare(df, options, common=None):
    """Apply the standard filters BEFORE selecting shared samples.

    Log-Pearson does not silently add a zero-row filter. Undefined correlations
    (constant/nonfinite rows, including invalid logs) cannot pass the standard
    correlation cutoff and can be excluded before block multiplication.
    """
    values = df.to_numpy(dtype=np.float64)
    keep = np.ones(len(df), dtype=bool)
    if options.mad_filter:
        keep &= np.median(np.abs(values - np.median(values, axis=1)[:, None]), axis=1) > 0
    if options.remove_zeros:
        keep &= (values != 0).all(axis=1)
    selected = df.loc[keep, common] if common is not None else df.loc[keep]
    values, names = selected.to_numpy(dtype=np.float64, copy=True), selected.index.to_numpy()
    if options.method == 'pearsonlog':
        with np.errstate(divide='ignore', invalid='ignore'):
            values = np.log10(values)
    elif options.method == 'spearman':
        values = np.asarray([stats.rankdata(row) for row in values])
    if not len(values):
        raise ValueError('No rows remain after filtering.')
    keep = np.isfinite(values).all(axis=1)
    values, names = values[keep], names[keep]
    values -= values.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(values, axis=1)
    keep = norms > 0
    log(f'Prepared {len(names)} rows; excluded {np.count_nonzero(~keep)} constant rows')
    if not keep.any():
        raise ValueError('All rows are constant across shared samples.')
    return np.ascontiguousarray(values[keep] / norms[keep, None]), names[keep]


def correlation_block(x, y, start, stop, cutoff, method, samples):
    r = np.clip(x[start:stop] @ y.T, -1.0, 1.0)
    i, j = np.nonzero(np.abs(r) >= cutoff)
    coeff = r[i, j]
    if samples == 2:
        p = np.full(len(coeff), np.nan) if method == 'spearman' else np.ones(len(coeff))
    elif method == 'spearman':
        with np.errstate(divide='ignore'):
            t = coeff * np.sqrt((samples - 2) / ((1 - coeff) * (1 + coeff)))
        p = 2 * stats.t.sf(np.abs(t), samples - 2)
    else:
        p = 2 * special.betainc(samples / 2 - 1, samples / 2 - 1, (1 - np.abs(coeff)) / 2)
    return start, stop, i + start, j, coeff, p


def blocks(x, y, options):
    """At most one in-flight/result block per worker; BLAS uses one thread each."""
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=options.threads) as executor:
        for base in range(0, len(x), options.block_size * options.threads):
            pending = []
            for start in range(base, min(len(x), base + options.block_size * options.threads), options.block_size):
                pending.append(executor.submit(correlation_block, x, y, start,
                                               min(start + options.block_size, len(x)),
                                               options.correlation_cutoff, options.method, x.shape[1]))
            for future in pending:
                yield future.result()


def append_frame(conn, table, frame):
    frame.to_sql(table, conn, if_exists='append', index=False, chunksize=10000)


def check_java(options):
    if not Path(options.clusterone_jar).is_file():
        raise ValueError(f'ClusterOne JAR missing: {options.clusterone_jar}. Set --clusterone-jar.')
    try:
        result = subprocess.run(['java', '-version'], capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("Java unavailable; activate/load Java and check 'java -version'.") from exc
    if result.returncode:
        raise RuntimeError("Java runtime unavailable: " + result.stderr)


def cluster(options, weights, outfile, output, table):
    try:
        with open(outfile, 'w') as handle:
            result = subprocess.run(['java', '-Xmx' + options.java_heap, '-jar', options.clusterone_jar,
                                     str(weights), '-F', 'csv'], stdout=handle, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError('Java unavailable when starting ClusterOne.') from exc
    if result.returncode:
        raise RuntimeError(f'ClusterOne exited {result.returncode}: {result.stderr.strip()}')
    if result.stderr.strip():
        log('ClusterOne: ' + result.stderr.strip())
    try:
        df = pd.read_csv(outfile)
    except pd.errors.EmptyDataError as exc:
        raise RuntimeError('ClusterOne produced empty output; inspect its diagnostics and network.') from exc
    append_frame(output, table, df)


def mutual_rank(scratch, output, options, nodes, directory):
    """Implement mutual_ranks.RANKING_CONVENTIONS with scratch SQLite.

    Node IDs are allocated in lexical partner order within each layer. Remove
    unimodal diagonals and validate/deduplicate pairs BEFORE ROW_NUMBER().
    """
    if sqlite3.sqlite_version_info < (3, 25, 0):
        raise RuntimeError('Mutual ranks require SQLite >= 3.25 (window functions).')
    with stage('Ranking retained edges on scratch disk'):
        unimodal = options.network_type == 'unimodal'
        left, right = ('min(a,b)', 'max(a,b)') if unimodal else ('a', 'b')
        where = 'WHERE a != b' if unimodal else ''
        scratch.execute(f'''CREATE TABLE cleaned AS
            SELECT {left} AS a, {right} AS b, min(correlation) AS low,
                   max(correlation) AS correlation, count(*) AS copies
            FROM edges {where} GROUP BY 1, 2''')
        conflicts = scratch.execute('''SELECT a,b,low,correlation FROM cleaned
            WHERE correlation-low > ? + ? * max(abs(low),abs(correlation)) LIMIT 10''',
            (mutual_ranks.DUPLICATE_ATOL, mutual_ranks.DUPLICATE_RTOL)).fetchall()
        if conflicts:
            raise ValueError(f'Conflicting duplicate correlation scores (first 10): {conflicts}')
        total = scratch.execute('SELECT count(*) FROM edges').fetchone()[0]
        self_rows = scratch.execute('SELECT count(*) FROM edges WHERE a=b').fetchone()[0] if unimodal else 0
        pairs = scratch.execute('SELECT count(*) FROM cleaned').fetchone()[0]
        options._rank_audit = dict(input_rows=total, self_rows_removed=self_rows,
                                  duplicate_rows_removed=total-self_rows-pairs, distinct_pairs=pairs)
        scratch.executescript('''
            CREATE TABLE directed AS
              SELECT a AS node, b AS neighbor, correlation FROM cleaned
              UNION ALL SELECT b, a, correlation FROM cleaned;
            DROP TABLE cleaned;
            DROP TABLE edges;
            CREATE INDEX directed_order ON directed(node, correlation DESC, neighbor);
            CREATE TABLE ranked AS
              SELECT node, neighbor, ROW_NUMBER() OVER
                (PARTITION BY node ORDER BY correlation DESC, neighbor) AS rank
              FROM directed;
            DROP TABLE directed;
            CREATE UNIQUE INDEX rank_pair ON ranked(node, neighbor);
            CREATE TABLE mr AS
              SELECT a.node AS a, a.neighbor AS b, a.rank * 1.0 * b.rank AS product
              FROM ranked a JOIN ranked b ON a.node=b.neighbor AND a.neighbor=b.node
              WHERE a.node < a.neighbor;
            DROP TABLE ranked;
        ''')
    write_mutual_ranks(pd.read_sql_query('SELECT a, b, product FROM mr ORDER BY a,b', scratch, chunksize=100000),
                       output, options, nodes, directory)


def write_mutual_ranks(chunks, output, options, nodes, directory):
    rates = options.multi_decay_rates or [options.decay_rate]
    empty = pd.DataFrame({'metabolite': pd.Series(dtype=str), 'gene': pd.Series(dtype=str),
                          'MR': pd.Series(dtype=float)})
    append_frame(output, options.sqlite_table_name + '_MR_edges', empty)
    for d in rates:
        suffix = f'_MR_weights_DR_{d}' if options.multi_decay_rates else '_MR_weights'
        append_frame(output, options.sqlite_table_name + suffix, empty.assign(weight=pd.Series(dtype=float)))
    files = {d: Path(directory) / f'weights_{d}.txt' for d in rates}
    handles = {d: open(path, 'w') for d, path in files.items()} if options.clusterone else {}
    counts = {d: 0 for d in rates}
    try:
        with stage('Writing mutual ranks and weights'):
            for chunk in chunks:
                if chunk.empty:
                    continue
                a, b = chunk.a.to_numpy(), chunk.b.to_numpy()
                mr = np.sqrt(chunk['product'].to_numpy())
                frame = pd.DataFrame({'metabolite': nodes[a], 'gene': nodes[b], 'MR': mr})
                append_frame(output, options.sqlite_table_name + '_MR_edges', frame)
                for d in rates:
                    weight = np.exp(-(mr - 1) / d)
                    keep = weight > options.edge_weight_cutoff
                    selected = frame.loc[keep].copy()
                    selected['weight'] = weight[keep]
                    suffix = f'_MR_weights_DR_{d}' if options.multi_decay_rates else '_MR_weights'
                    append_frame(output, options.sqlite_table_name + suffix, selected)
                    counts[d] += int(keep.sum())
                    if options.clusterone:
                        mutual_ranks.clusterone_edges(selected, network_type=options.network_type).to_csv(
                            handles[d], sep=' ', header=False, index=False)
    finally:
        for handle in handles.values():
            handle.close()
    if options.clusterone:
        for d in rates:
            if not counts[d]:
                log(f'No retained edges at decay {d}; skipping ClusterOne')
                continue
            with stage(f'ClusterOne decay={d} edges={counts[d]}'):
                suffix = f'_clone_DR_{d}' if options.multi_decay_rates else '_clone'
                cluster(options, files[d], Path(directory) / f'clusters_{d}.csv',
                        output, options.sqlite_table_name + suffix)



EDGE_DTYPE = np.dtype([('a', '<i8'), ('b', '<i8'), ('correlation', '<f8')])
RANK_DTYPE = np.dtype([('node', '<i8'), ('neighbor', '<i8'), ('score', '<f8'), ('rank', '<i8')])


def ram_rank_estimate(edge_count):
    # Conservative allowance for input chunks, directed arrays, deduplication,
    # sorting/selection temporaries, and fixed output buffers. Not an RSS limit.
    return edge_count * 256 + 128 * 2**20


def check_ram_budget(edge_count, budget_gb):
    estimated = ram_rank_estimate(edge_count) / 2**30
    if budget_gb is not None and estimated > budget_gb:
        raise ValueError(f'RAM mutual ranks need an estimated {estimated:.2f} GiB for '
                         f'{edge_count:,} retained edges, exceeding --mr-memory-gb {budget_gb:.2f}. '
                         'Use --mr-storage disk or increase the budget and SLURM memory allocation. '
                         'This run may have partial output tables; use a fresh -tn when rerunning.')
    return estimated


def rank_edges_ram(edge_chunks, network_type):
    """Rank numeric keys and scatter ranks back to canonical edge order.

    lexsort works on primitive numeric columns; sorting structured records by
    several fields repeatedly is substantially slower at millions of edges.
    Keep each reciprocal edge at its original array position so neither a
    second rank-order sort nor a reverse-edge join is needed.
    """
    count = sum(len(chunk) for chunk in edge_chunks)
    edges = np.empty(count, dtype=EDGE_DTYPE)
    offset = 0
    for chunk in edge_chunks:
        stop = offset + len(chunk)
        edges[offset:stop] = chunk
        offset = stop
    edge_chunks.clear()
    if network_type == 'unimodal':
        edges = edges[edges['a'] != edges['b']]
        a = np.minimum(edges['a'], edges['b'])
        edges['b'] = np.maximum(edges['a'], edges['b'])
        edges['a'] = a
        del a
    with stage(f'RAM pair cleanup: {len(edges):,} non-self input edges'):
        order = np.lexsort((edges['b'], edges['a']))
        edges = edges[order]
        del order
        if len(edges):
            starts = np.r_[0, np.flatnonzero((edges['a'][1:] != edges['a'][:-1]) |
                                            (edges['b'][1:] != edges['b'][:-1])) + 1]
            low = np.minimum.reduceat(edges['correlation'], starts)
            high = np.maximum.reduceat(edges['correlation'], starts)
            conflicts = high-low > mutual_ranks.DUPLICATE_ATOL + mutual_ranks.DUPLICATE_RTOL * np.maximum(abs(low), abs(high))
            if conflicts.any():
                examples = [(int(edges['a'][i]), int(edges['b'][i]), float(lo), float(hi))
                            for i, lo, hi in zip(starts[conflicts][:10], low[conflicts][:10], high[conflicts][:10])]
                raise ValueError(f'Conflicting duplicate correlation scores (first 10): {examples}')
            edges = edges[starts]
            edges['correlation'] = high  # same deterministic representative as the disk backend
            del starts, low, high, conflicts
    count = len(edges)
    directed = np.empty(2 * count, dtype=RANK_DTYPE)
    directed['node'][:count], directed['neighbor'][:count] = edges['a'], edges['b']
    directed['node'][count:], directed['neighbor'][count:] = edges['b'], edges['a']
    directed['score'][:count] = directed['score'][count:] = -edges['correlation']
    directed['rank'] = 0
    del edges
    with stage(f'RAM partner ranking: {count:,} distinct pairs'):
        if count:
            order = np.lexsort((directed['neighbor'], directed['score'], directed['node']))
            ordered_nodes = directed['node'][order]
            starts = np.r_[0, np.flatnonzero(ordered_nodes[1:] != ordered_nodes[:-1]) + 1]
            sizes = np.diff(np.r_[starts, len(order)])
            ranks = np.arange(1, len(order)+1, dtype=np.int64) - np.repeat(starts, sizes)
            directed['rank'][order] = ranks
            del order, ordered_nodes, starts, sizes, ranks
    # Both halves still follow the same canonical pair order. Restore canonical
    # node/neighbor columns in the reverse half, keeping the reverse ranks.
    swap = directed['node'][count:].copy()
    directed['node'][count:] = directed['neighbor'][count:]
    directed['neighbor'][count:] = swap
    return directed[:count], directed[count:]


def ram_rank_chunks(edge_chunks, options):
    # Own arrays inside this generator: exhaustion releases them before Java starts.
    count = sum(len(chunk) for chunk in edge_chunks)
    estimated = check_ram_budget(count, options.mr_memory_gb)
    log(f'RAM mutual ranks: retained edges={count:,}; estimated working memory={estimated:.2f} GiB '
        '(additional process memory is required)')
    with stage('Ranking retained edges in RAM'):
        self_rows = sum(int((chunk['a'] == chunk['b']).sum()) for chunk in edge_chunks) if options.network_type == 'unimodal' else 0
        forward, backward = rank_edges_ram(edge_chunks, options.network_type)
        options._rank_audit = dict(input_rows=count, self_rows_removed=self_rows,
                                  duplicate_rows_removed=count-self_rows-len(forward), distinct_pairs=len(forward))
    for start in range(0, len(forward), 100000):
        a, b = forward[start:start+100000], backward[start:start+100000]
        yield pd.DataFrame({'a': a['node'], 'b': a['neighbor'],
                            'product': a['rank'].astype(np.float64) * b['rank']}, copy=True)


def ram_mutual_rank(edge_chunks, output, options, nodes, directory):
    write_mutual_ranks(ram_rank_chunks(edge_chunks, options), output, options, nodes, directory)


def main(options):
    # Match the standard runner's protection of all original MR databases.
    if options.mutual_rank and Path(options.sqlite_db_name).exists():
        raise FileExistsError('MR construction requires a new output database; preserve the existing database.')
    options.network_type = mutual_ranks.network_type_from_inputs(options.feature_table, options.quantitation_matrix)
    do_mr = options.mutual_rank and options.method != 'spearman'
    if options.mr_storage not in ('ram', 'disk'):
        raise ValueError('--mr-storage must be ram or disk')
    log(f'Runner={Path(__file__).resolve()}; SHA256={file_sha256(__file__)}')
    log(f'Network type={options.network_type}; feature input={Path(options.feature_table).resolve()}; '
        f'transcript input={Path(options.quantitation_matrix).resolve()}')
    if options.mr_memory_gb is not None and (not np.isfinite(options.mr_memory_gb) or options.mr_memory_gb <= 0):
        raise ValueError('--mr-memory-gb must be a positive finite number (GiB).')
    if options.mr_storage == 'ram' and options.mr_memory_gb is None:
        allocated_mb = os.environ.get('SLURM_MEM_PER_NODE')
        if allocated_mb and allocated_mb.isdigit() and int(allocated_mb) > 0:
            options.mr_memory_gb = 0.6 * int(allocated_mb) / 1024
    if options.mutual_rank:
        log(f'Mutual-rank storage={options.mr_storage}; RAM working budget='
            f'{options.mr_memory_gb if options.mr_memory_gb is not None else "not set"} GiB')
    allocation = int(os.environ.get('SLURM_CPUS_PER_TASK', options.threads))
    if options.threads < 1 or options.threads > allocation or options.block_size < 1:
        raise ValueError('Use positive threads/block-size; threads must not exceed SLURM_CPUS_PER_TASK.')
    if not 0 <= options.correlation_cutoff <= 1 or not 0 <= options.edge_weight_cutoff <= 1:
        raise ValueError('Correlation/weight cutoffs must be between 0 and 1.')
    rates = options.multi_decay_rates or [options.decay_rate]
    if any(d <= 0 for d in rates) or len(set(rates)) != len(rates):
        raise ValueError('Decay rates must be positive and unique.')
    if options.clusterone and do_mr:
        check_java(options)
    with stage('Read and prepare inputs'):
        transcripts, features = read_matrix(options.quantitation_matrix), read_matrix(options.feature_table)
        common = sorted(set(transcripts.columns) & set(features.columns))
        if len(common) < 2:
            raise ValueError(f'Need >=2 shared samples; found {len(common)}. Transcript headers: '
                             f'{list(transcripts.columns[:5])}; feature headers: {list(features.columns[:5])}')
        if len(common) < min(len(transcripts.columns), len(features.columns)):
            log(f'WARNING: Using {len(common)} shared samples; unmatched columns excluded')
        x, feature_names = prepare(features, options, common)
        y, gene_names = prepare(transcripts, options, common)
        if options.network_type == 'unimodal':
            nodes = np.asarray(sorted(set(feature_names) | set(gene_names)), dtype=object)
            lookup = {node: i for i, node in enumerate(nodes)}
            feature_ids = np.asarray([lookup[n] for n in feature_names])
            gene_ids = np.asarray([lookup[n] for n in gene_names])
            layers = np.repeat('feature', len(nodes))
            options.cluster_nodes = nodes
        else:
            feature_layer, gene_layer = sorted(feature_names), sorted(gene_names)
            nodes = np.asarray(feature_layer + gene_layer, dtype=object)
            layers = np.asarray(['metabolite'] * len(feature_layer) + ['gene'] * len(gene_layer))
            feature_lookup = {name: i for i, name in enumerate(feature_layer)}
            gene_lookup = {name: len(feature_layer) + i for i, name in enumerate(gene_layer)}
            feature_ids = np.asarray([feature_lookup[n] for n in feature_names])
            gene_ids = np.asarray([gene_lookup[n] for n in gene_names])
            options.cluster_nodes = np.asarray([layer + ':' + escape_identifier(name, safe='')
                                                for layer, name in zip(layers, nodes)], dtype=object)
        log(f'features={len(x)} transcripts={len(y)} samples={len(common)} pairs={len(x)*len(y):,} '
            f'workers={options.threads} block={options.block_size}')
        log(f'Normalized matrices={((x.nbytes+y.nbytes)/2**20):.1f} MiB; '
            f'raw correlation blocks <= {options.threads*options.block_size*len(y)*8/2**20:.1f} MiB '
            '(additional index, p-value and SQL buffers required)')
        del transcripts, features
    # Job-unique scratch files and a single database writer avoid concurrent worker writes.
    with tempfile.TemporaryDirectory(prefix='corrMO_', dir=options.scratch_dir) as directory, ExitStack() as connections:
        rank_path = Path(directory) / 'ranks.sqlite'
        scratch = None
        if do_mr and options.mr_storage == 'disk':
            scratch = sqlite3.connect(str(rank_path))
            connections.callback(scratch.close)
            connections.enter_context(scratch)
        output = sqlite3.connect(options.sqlite_db_name, timeout=60)
        connections.callback(output.close)
        connections.enter_context(output)
        if do_mr:
            if options.mr_storage == 'disk':
                log(f'DISK ranking database={rank_path}')
            else:
                log('RAM ranking uses NumPy arrays; no scratch ranking database is opened.')
            log(f'Persistent correlation/MR/weight output={Path(options.sqlite_db_name).resolve()}; '
                f'ClusterONE intermediates={directory}')
        if scratch is not None:
            scratch.execute('PRAGMA cache_size=-262144')
            scratch.execute('PRAGMA temp_store=FILE')
            # Only disposable, job-local ranks use these settings. Leave the
            # persistent output database's durability settings untouched.
            scratch.execute('PRAGMA journal_mode=OFF')
            scratch.execute('PRAGMA synchronous=OFF')
            scratch.execute('CREATE TABLE edges (a INTEGER, b INTEGER, correlation REAL)')
        edge_chunks = []
        prefix = options.sqlite_table_name
        tables = [prefix, prefix + '_nodes']
        if do_mr:
            tables += [prefix + '_MR_edges', prefix + '_MR_provenance']
            tables += [prefix + f'_MR_weights_DR_{d}' for d in rates] if options.multi_decay_rates else [prefix + '_MR_weights']
            if options.clusterone:
                tables += [prefix + f'_clone_DR_{d}' for d in rates] if options.multi_decay_rates else [prefix + '_clone']
        existing = {r[0] for r in output.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        collisions = existing.intersection(tables)
        if collisions and options.if_exists == 'fail':
            raise ValueError(f'Output tables already exist: {sorted(collisions)}. Use a new -tn, '
                             '--if-exists replace, or explicitly opt into append.')
        if options.if_exists == 'replace':
            for table in collisions:
                output.execute('DROP TABLE ' + quote(table))
            output.commit()
        mapping = pd.DataFrame({'node_id': np.arange(len(nodes)), 'label': nodes, 'layer': layers,
                                'clusterone_label': options.cluster_nodes})
        if prefix + '_nodes' in collisions and options.if_exists == 'append':
            previous = pd.read_sql_query('SELECT * FROM ' + quote(prefix + '_nodes'), output)
            if not previous.equals(mapping):
                raise ValueError('Cannot append: node ID mapping differs from the existing run. Use a new -tn.')
        else:
            append_frame(output, prefix + '_nodes', mapping)
        retained = 0
        append_frame(output, prefix, pd.DataFrame({'metabolite': pd.Series(dtype=str),
            'gene': pd.Series(dtype=str), 'correlation': pd.Series(dtype=float), 'P': pd.Series(dtype=float)}))
        with stage('Block correlations and incremental database writes'):
            for start, stop, i, j, r, p in blocks(x, y, options):
                append_frame(output, prefix, pd.DataFrame({'metabolite': feature_names[i],
                             'gene': gene_names[j], 'correlation': r, 'P': p}))
                if do_mr:
                    if options.mr_storage == 'disk':
                        scratch.executemany('INSERT INTO edges VALUES (?, ?, ?)',
                                            zip(feature_ids[i].tolist(), gene_ids[j].tolist(), r.tolist()))
                        scratch.commit()
                    else:
                        check_ram_budget(retained + len(r), options.mr_memory_gb)
                        edge = np.empty(len(r), dtype=EDGE_DTYPE)
                        edge['a'], edge['b'], edge['correlation'] = feature_ids[i], gene_ids[j], r
                        edge_chunks.append(edge)
                retained += len(r)
                log(f'Correlation rows {stop}/{len(x)}; retained pairs={retained:,}')
        del x, y  # release normalized matrices before allocating ranking arrays
        if do_mr:
            with stage('Mutual rank pipeline'):
                if options.mr_storage == 'ram':
                    ram_mutual_rank(edge_chunks, output, options, nodes, directory)
                else:
                    mutual_rank(scratch, output, options, nodes, directory)
            settings = {key: value for key, value in vars(options).items()
                        if key not in ('cluster_nodes', '_rank_audit')}
            # Programmatic callers may use Path arguments rather than CLI strings.
            settings = {key: str(value) if isinstance(value, Path) else value for key, value in settings.items()}
            provenance = mutual_ranks.ranking_provenance(options.network_type, options=settings)
            provenance.update(input_audit=options._rank_audit,
                network_type_selection='same input file: unimodal; different input files: bipartite',
                caller_sha256=file_sha256(__file__),
                input_sha256={str(Path(path).resolve()): file_sha256(path)
                              for path in (options.feature_table, options.quantitation_matrix)})
            if options.clusterone:
                provenance['clusterone_jar_sha256'] = file_sha256(options.clusterone_jar)
                provenance['clusterone_arguments'] = ['-Xmx' + options.java_heap, '-F', 'csv']
            append_frame(output, prefix + '_MR_provenance', pd.DataFrame([{
                'created_utc': datetime.now(timezone.utc).isoformat(),
                'provenance_json': json.dumps(provenance, sort_keys=True)}]))
        elif not retained:
            log('No correlations passed cutoff; mutual ranks and clustering skipped')
    log('Completed successfully')


if __name__ == '__main__':
    try:
        main(get_args())
    except MemoryError:
        print('ERROR: Insufficient RAM. Retry with --mr-storage disk or request more SLURM memory. '
              'Use a fresh -tn because output tables may be partial.', file=sys.stderr, flush=True)
        sys.exit(1)
    except (ValueError, RuntimeError, OSError, sqlite3.Error) as exc:
        print(f'ERROR: {exc}', file=sys.stderr, flush=True)
        sys.exit(1)
