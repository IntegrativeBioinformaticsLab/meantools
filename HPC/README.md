# HPC correlation and mutual-rank runner

`corrMultiomics_HPC.py` follows the correlation, filtering, ranking, weighting,
and ClusterONE conventions of the repository's `corrMultiomics.py`. Its
implementation uses block correlations, a bounded thread pool, streamed SQLite
writes, and RAM or scratch-disk ranks for single-node SLURM jobs.

Keep the `HPC/` folder inside the full MEANtools checkout. The runner imports the
shared `mutual_ranks.py` from the parent directory for automatic layer selection,
ranking conventions, duplicate tolerances, ClusterONE labels and provenance.
It does not import `gizmos` or require rdkit for correlation/network construction.

## Shared analysis behavior

- Use the same file for `-ft` and `-qm` for gene–gene coexpression or
  metabolite–metabolite coabundance. Different files define a gene–metabolite
  bipartite network, even when identifiers or complete matrix contents match.
  No network-type parameter is exposed. Symlinks and hard links to the same
  input are recognized as unimodal.
- MAD and optional zero-row filters use the original input columns, before
  selecting the sorted intersection of sample columns, matching the standard
  runner. Numerical centering and scaling implement Pearson correlation; no
  biological input normalization is added.
- `pearsonlog` applies `log10` after filtering/sample selection. It does not add
  a separate zero-row filter. Rows with undefined correlations, including
  constants, invalid logarithms and nonfinite values, cannot pass the standard
  correlation cutoff and are omitted from block multiplication.
- Pearson, log-Pearson and Spearman coefficients/p-values match SciPy within
  floating-point precision. As in the standard runner, Spearman writes
  correlations but skips MR/weights/ClusterONE even when `-mr` is supplied.
- Correlation filtering remains `abs(correlation) >= cutoff`. Partner ranking
  uses signed correlation descending. Exact ties use partner identifier
  ascending, with ordinal ranks starting at 1. Internal numeric IDs are
  allocated in lexical order within each biological layer.
- Both ranking backends remove unimodal self-correlations and canonicalize
  undirected pairs before assigning ranks. Bipartite nodes retain separate
  layer identities and equal-name gene–metabolite edges remain eligible.
- Duplicate score extrema must satisfy
  `max-min <= 1e-12 + 1e-9 * max(abs(min),abs(max))`. Conflicts raise an error;
  agreeing duplicates use the maximum score, matching `mutual_ranks.py`.
- MR is `sqrt(rank_A(B) * rank_B(A))`. Weights remain
  `exp(-(MR - 1) / decay_rate)` with the same strict `weight > cutoff` rule.
  Default decay rates remain `[5, 10, 25]` and the edge cutoff remains `0.01`.
- ClusterONE receives the same ordered labels, edges, weights and clustering
  arguments as the standard runner. Bipartite labels are URL-escaped and
  prefixed by `metabolite:` / `gene:`. Module tables use the standard schema;
  numeric scratch IDs are not substituted for biological labels at clustering.
  The HPC Java heap argument controls resources, not clustering parameters.
- `-p` and `-a` retain the standard runner's current no-op behavior.
- MR construction requires a **new output database**, including when the HPC
  `--if-exists` argument is supplied. Original databases are preserved.

Block matrix products can differ from SciPy's pairwise arithmetic at the last
floating-point digits. Exact numerical ties and scores exactly on a cutoff
boundary can therefore be sensitive to that arithmetic. Both rank backends
apply the documented exact-tie policy to their actual retained scores.

## HPC execution differences

The runner needs Python >=3.8, NumPy, SciPy, pandas, `threadpoolctl`, and SQLite
>=3.25 for disk ranks. ClusterONE additionally needs Java and the repository's
`cluster_one-1.0.jar`. Install `threadpoolctl` in the HPC environment if absent.

| Option | Purpose |
| --- | --- |
| `-t` | Correlation workers; defaults to `SLURM_CPUS_PER_TASK`, or 1 outside SLURM |
| `--block-size` | Feature rows per block; default 128 |
| `--mr-storage disk` | Default: group, validate and rank retained edges in scratch SQLite |
| `--mr-storage ram` | Group, validate and rank compact NumPy edge arrays in memory |
| `--mr-memory-gb` | RAM-mode estimated working-memory budget in GiB |
| `--scratch-dir` | Job-local intermediates; defaults to `SLURM_TMPDIR` when available |
| `--clusterone-jar` | JAR location; discovers the repository's bundled JAR by default |
| `--java-heap` | Java maximum heap; default `8g` |
| `--if-exists` | Table collision policy for correlation-only runs; MR runs always require a new database |

One BLAS thread per worker avoids nested CPU oversubscription. At most one
correlation block per worker is held in flight. There are no multiprocessing
managers, forked worker copies, or concurrent database writers.

Both MR backends rank retained edges directly, avoiding all-node combinations.
The disk backend uses SQL grouping, window ranks and a reciprocal join. The RAM
backend sorts numeric column keys, deduplicates before rank assignment and
scatters partner ranks back to canonical edge positions. This avoids repeated
structured-record sorts and a separate reverse-edge join. RAM planning uses
the existing estimate of 256 bytes per
retained input edge plus 128 MiB; this is an estimate, not a hard RSS limit.
Input edge counts include reverse and diagonal records before cleanup. If no
explicit RAM budget is set, a positive numeric `SLURM_MEM_PER_NODE` supplies a
budget of 60% of the allocation. Ranking arrays are released before Java starts.

SQLite intermediate storage, cache limits, block size, worker limits and Java
heap settings affect resource use. Scratch sorting may need substantial local
space. Stage timings and block progress are logged. The reader additionally
recognizes semicolon text CSVs and reports malformed/archived input errors.
No GPU or multi-node distribution is implemented.

## Selecting RAM mode and diagnosing disk activity

The submission wrapper passes the storage mode exactly once, using
`MR_STORAGE` (default `disk`). Select RAM when submitting:

```bash
sbatch --export=ALL,MR_STORAGE=ram HPC/corrMultiomics.sh
```

For a direct Python invocation, supply `--mr-storage ram` once. In older
wrappers, a command containing `--mr-storage ram` followed by `"${mr_args[@]}"`
could also include the array's default `--mr-storage disk`. Python's default
argument handling silently used that later value. The updated parser rejects
conflicting repeated storage options before reading data or creating output.
Update both HPC files in the cluster installation before resubmitting.

Startup logs identify the actual runner path and SHA256, selected network type,
effective storage mode, memory budget and output paths. In RAM mode, the log
states that NumPy arrays are used and no scratch ranking database is opened.
Only disk mode opens `ranks.sqlite`. A growing `ranks.sqlite` therefore indicates
disk ranking or an older installed script, rather than RAM automatically
spilling to disk. Exceeding the RAM budget raises an error; it does not change
backends.

RAM mode still writes persistent correlations, MRs and weights to the output
SQLite database and creates edge files for ClusterONE. Those output writes
remain disk I/O. For disk ranking, use node-local scratch when your cluster
provides it; if a path under `/nobackup/beegfs` is shared storage, large SQLite
sorts and intermediate writes can be limited by that filesystem. Journaling
and synchronous writes are disabled only for the disposable scratch ranking
database; output-database durability settings are preserved.

For unimodal coexpression, pass the **same chosen transcript file** to both
`-ft` and `-qm`. Two separately named files, such as `expression_2.csv` and
`expression.csv`, select bipartite handling even if their contents are equal.
The startup network-type log makes this distinction visible.

A local ranking-only benchmark with 3,000 nodes and 7,201,814 synthetic input
records took 18.6 seconds with the previous RAM sorter and 5.3 seconds with the
numeric-key implementation (about 3.5-fold faster). Forward and reciprocal
rank arrays were exactly equal. These timings exclude correlation calculation,
SQLite output and ClusterONE and do not predict full-job runtime on SLURM.

## Outputs

The standard table names and schemas are retained:

- `<prefix>`: `metabolite, gene, correlation, P`.
- `<prefix>_MR_edges`: `metabolite, gene, MR`.
- `<prefix>_MR_weights_DR_<d>` and `<prefix>_clone_DR_<d>` for multi-decay runs.
- `<prefix>_MR_weights` and `<prefix>_clone` for the programmatic single-decay path.
- `<prefix>_MR_provenance`: timestamp, conventions, automatic layer selection,
  cleanup counts, software versions, input/caller hashes and JAR settings.

The additional HPC `<prefix>_nodes` table records internal numeric IDs, raw
labels, layers and ClusterONE labels. It is metadata for scratch representation;
MR/weight tables keep original endpoint identifiers. Empty MR and weight tables
are still created when no non-self pair remains. ClusterONE is skipped when no
edge survives its weight cutoff. Writes are incremental: a failed job can leave
a partial output database; restart with a new database path.

## Submit

Review the site-specific paths, modules, account and resource directives in
`HPC/corrMultiomics.sh` before submission:

```bash
sbatch HPC/corrMultiomics.sh
sbatch --export=ALL,MR_STORAGE=ram,MR_MEMORY_GB=20 HPC/corrMultiomics.sh
```

The example now uses the same transcript matrix for both inputs, so its
coexpression network is unimodal, and includes the SLURM job ID in the output
**database filename**. For bipartite analysis, change `-ft` to the separate
metabolite matrix. Intermediates use job-local scratch and are cleaned on exit.

For direct use from any working directory:

```bash
python /path/to/meantools/HPC/corrMultiomics_HPC.py \
  -ft /path/to/transcripts.csv -qm /path/to/transcripts.csv \
  -dn /path/to/new_coexpression.sqlite -tn coexpression \
  -mr -cl -c 0.1 -w 0.01 -mdr 5 10 25 \
  -t 8 --mr-storage disk --scratch-dir /path/to/scratch
```

Use the same input matrices and analysis options when comparing standard and
HPC runs. Inspect SLURM accounting and stage logs to assess full-job resources;
local tests do not establish performance or queue behavior on the cluster.

## Validation

```bash
python -m unittest discover -s HPC/tests -v
python -m pytest tests HPC/tests -q
bash -n HPC/corrMultiomics.sh
```

Tests compare both rank backends with `mutual_ranks.py` for upper triangles,
symmetric/repeated/shuffled edges, self-pairs, tied/negative correlations,
duplicate tolerance/conflicts, empty graphs and overlapping bipartite names.
They also exercise filtering order, log handling, streaming, provenance,
RAM budgets, database preservation and Spearman's correlation-only behavior.
Storage-selection tests execute the wrapper with site commands stubbed, reject
conflicting options, and trace database connections to verify that RAM mode
opens only the persistent output database and never the scratch rank database.
In a full MEANtools environment, integration tests run the actual standard and
HPC CLIs, compare correlations/MRs/weights, and compare real ClusterONE member
sets for both network types, both storage backends and both decay-rate paths.
The integration test is skipped if rdkit is unavailable for the standard runner.
