#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import hashlib
import platform
from pathlib import Path
import subprocess as sp
from urllib.parse import quote
import pandas as pd
import numpy as np

DUPLICATE_ATOL = 1e-12
DUPLICATE_RTOL = 1e-9
RANKING_CONVENTIONS = {
    'version': 'distinct-partners-v1',
    'order': 'signed correlation descending (not absolute correlation)',
    'ties': 'ordinal ranks starting at 1; exact ties by partner identifier ascending',
    'duplicates': 'validate extrema: max-min <= atol + rtol*max(abs(min),abs(max)); retain max',
    'duplicate_atol': DUPLICATE_ATOL,
    'duplicate_rtol': DUPLICATE_RTOL,
    'unimodal': 'canonical undirected pairs; exclude self before ranking',
    'bipartite': 'column defines node layer; equal identifiers across layers are distinct',
    'MR': 'sqrt(rank_A(B) * rank_B(A))',
}



def calc_MR(rank1, rank2):
    MR = np.sqrt(rank1 * rank2)
    return MR


def network_type_from_inputs(feature_table, quantitation_matrix):
    """Use existing input-file identity to select node namespaces automatically.

    Comparing one matrix with itself defines a unimodal network. Different
    input files define separate layers, even if their identifiers or contents
    match. samefile also recognizes relative paths, symlinks and hard links.
    Correlation rows alone cannot distinguish these biological namespaces.
    """
    return ('unimodal' if Path(feature_table).samefile(quantitation_matrix)
            else 'bipartite')


def clean_correlations(correlations_df, *, network_type, atol=DUPLICATE_ATOL,
                       rtol=DUPLICATE_RTOL):
    """Return one validated score per pair, without changing correlation filtering.

    Identifiers must be strings. Unimodal endpoints share a namespace and are
    canonicalized lexicographically. Bipartite endpoints are layer-specific:
    metabolite A / gene B and metabolite B / gene A are different edges.
    Duplicate score extrema must agree within the symmetric tolerance documented
    in RANKING_CONVENTIONS. The maximum is a deterministic representative, so
    repeated rows do not change it. Tolerance is NOT used to define rank ties.
    """
    if network_type not in ('unimodal', 'bipartite'):
        raise ValueError("network_type must be 'unimodal' or 'bipartite'")
    if not np.isfinite([atol, rtol]).all() or atol < 0 or rtol < 0:
        raise ValueError('Duplicate tolerances must be finite and nonnegative')
    df = correlations_df[['metabolite', 'gene', 'correlation']].copy()
    for column in ('metabolite', 'gene'):
        if not df[column].map(lambda value: isinstance(value, str)).all():
            raise ValueError('Node identifiers must be non-null strings')
    df['correlation'] = pd.to_numeric(df.correlation, errors='raise')
    if not np.isfinite(df.correlation.to_numpy(dtype=float)).all():
        raise ValueError('Correlation scores must be finite')
    input_rows = len(df)
    self_rows = 0
    if network_type == 'unimodal':
        self_mask = df.metabolite == df.gene
        self_rows = int(self_mask.sum())
        df = df.loc[~self_mask].copy()
        swap = df.metabolite > df.gene
        df.loc[swap, ['metabolite', 'gene']] = df.loc[swap, ['gene', 'metabolite']].to_numpy()
    grouped = df.groupby(['metabolite', 'gene'], sort=True).correlation.agg(['min', 'max', 'size'])
    tolerance = atol + rtol * np.maximum(grouped['min'].abs(), grouped['max'].abs())
    conflicts = grouped[grouped['max'] - grouped['min'] > tolerance]
    if len(conflicts):
        raise ValueError(
            f'Conflicting duplicate correlation scores for {len(conflicts)} pair(s) '
            f'(atol={atol}, rtol={rtol}); first 10:\n{conflicts.head(10).to_string()}'
        )
    cleaned = grouped['max'].rename('correlation').reset_index()
    cleaned.attrs['audit'] = dict(input_rows=input_rows, self_rows_removed=self_rows,
                                  duplicate_rows_removed=len(df) - len(cleaned),
                                  distinct_pairs=len(cleaned),
                                  pairs_with_duplicate_rows=int((grouped['size'] > 1).sum()),
                                  min_pair_multiplicity=int(grouped['size'].min()) if len(grouped) else 0,
                                  max_pair_multiplicity=int(grouped['size'].max()) if len(grouped) else 0,
                                  max_duplicate_score_difference=float((grouped['max'] - grouped['min']).max()) if len(grouped) else 0.)
    return cleaned


def get_MR_from_correlations(correlations_df, all_objects=None, *, network_type,
                             atol=DUPLICATE_ATOL, rtol=DUPLICATE_RTOL):
    """Recompute reciprocal ranks from distinct partners (see RANKING_CONVENTIONS).

    Network type is required; never infer biological layers from text overlap.
    The optional legacy all_objects argument must cover all input endpoints;
    it does not determine bipartite node identity. Output retains the historical
    column names, with canonical endpoints for unimodal networks and original
    layers for bipartite networks. Input data is never modified.
    """
    df = clean_correlations(correlations_df, network_type=network_type, atol=atol, rtol=rtol)
    if all_objects is not None:
        missing = (set(df.metabolite) | set(df.gene)) - set(all_objects)
        if missing:
            raise ValueError(f'all_objects is missing input endpoints: {sorted(missing)[:10]}')
    # Integer edge indices align reciprocal ranks without an all-node O(V^2)
    # combinations set. Layer identity is preserved by ranking each side alone.
    def rank_side(edges, node, partner):
        ordered = edges.sort_values([node, 'correlation', partner], ascending=[True, False, True])
        return (ordered.groupby(node, sort=False).cumcount() + 1).reindex(df.index)

    if network_type == 'bipartite':
        rank1 = rank_side(df, 'metabolite', 'gene')
        rank2 = rank_side(df, 'gene', 'metabolite')
    else:
        forward = df.rename(columns={'metabolite': 'node', 'gene': 'partner'})
        reverse = df.rename(columns={'gene': 'node', 'metabolite': 'partner'})
        incident = pd.concat([forward, reverse], ignore_index=True)
        ordered = incident.sort_values(['node', 'correlation', 'partner'], ascending=[True, False, True])
        ranks = (ordered.groupby('node', sort=False).cumcount() + 1).reindex(incident.index)
        rank1 = ranks.iloc[:len(df)].to_numpy(dtype=float)
        rank2 = ranks.iloc[len(df):].to_numpy(dtype=float)
    edges = df[['metabolite', 'gene']].copy()
    edges['MR'] = calc_MR(rank1, rank2)
    edges.attrs['audit'] = df.attrs['audit']
    return edges


def clusterone_edges(weights_df, *, network_type):
    """Export layer-safe ClusterONE nodes; bipartite tokens are URL-escaped.

    Raw identifiers remain in MR/weight tables. All bipartite module members use
    metabolite:<escaped-id> / gene:<escaped-id>, even when the names are disjoint.
    """
    if network_type not in ('unimodal', 'bipartite'):
        raise ValueError('Invalid network_type')
    edges = weights_df[['metabolite', 'gene', 'weight']].copy()
    if network_type == 'bipartite':
        for layer in ('metabolite', 'gene'):
            edges[layer] = edges[layer].map(lambda name: layer + ':' + quote(name, safe=''))
    return edges


def ranking_provenance(network_type, **settings):
    """Serializable conventions plus exact implementation and environment."""
    source = Path(__file__).resolve()
    git = sp.run(['git', '-C', str(source.parent), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    return dict(ranking=RANKING_CONVENTIONS, network_type=network_type,
                implementation_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                git_head=git.stdout.strip() if git.returncode == 0 else None,
                python=platform.python_version(), pandas=pd.__version__, numpy=np.__version__,
                settings=settings)

def get_weight_from_MR(edges_df, edge_weight_cutoff, decay_rate):
    # CONVERT TO EDGE WEIGHT
    print('Estimating edge weights using exponetial decay...')
    edges_df['weight'] = np.exp(-(edges_df.MR - 1) / decay_rate)  # other options: /5 /10 /25 /50 /100
    weights_df = edges_df[edges_df.weight > edge_weight_cutoff].reset_index(drop=True) # default edge_weight_cutoff is 0.01
    
    return weights_df

def run_clusterone(outfile, clusterone, weights_file):
    # CLUSTERONE
    print('\n')
    with open(outfile, 'w') as f:
        # java -jar cluster_one-1.0.jar input -F csv
        sp.run(['java', '-jar', clusterone, weights_file, '-F', 'csv'], stdout=f, check=True)
    print('\n')

    return

