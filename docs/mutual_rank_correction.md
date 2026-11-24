# Mutual-rank correction and network comparison

The ranking bug is fixed in `mutual_ranks.py`. Unimodal reverse records and
self-correlations previously consumed extra rank positions; a later dictionary
assignment retained an inflated rank for each repeated partner. The correction
recomputes ranks from correlations, rather than rescaling stored MRs.

## Ranking conventions

- `corrMultiomics.py` automatically selects a unimodal network when the
  existing `-ft` and `-qm` paths refer to the same file, including symlinks and
  hard links. Different input files define a bipartite network. There is no
  user-facing network-type option. The caller supplies the selected type to
  the ranking/export helpers internally and records the choice in provenance.
  For a unimodal run, use the same matrix path for both inputs rather than
  making a separate copy. Different files remain separate layers even when
  feature identifiers or complete matrix contents match; neither text overlap
  nor numeric equality establishes biological node identity.
- In a unimodal network, discard diagonal records, then canonicalize each pair
  by lexicographically sorting its two string identifiers.
- Validate the minimum and maximum duplicate score using
  `max - min <= 1e-12 + 1e-9 * max(abs(min), abs(max))`. Raise `ValueError` on a
  conflict, including the number of conflicting pairs and the first ten pairs'
  score ranges. For agreeing duplicates, use the maximum score. This choice is
  independent of record order and repeated copies; do not average duplicates.
- Rank distinct partners by **signed correlation descending**, as before.
  Negative correlations remain eligible under the existing absolute-correlation
  filter and sort below positive correlations. No filter is added here.
- Use ordinal integer ranks starting at 1. Break exact correlation ties by
  partner identifier ascending. Duplicate tolerance does not define a tie.
- Compute `MR(A,B) = sqrt(rank_A(B) * rank_B(A))` once per retained pair.
- Bipartite endpoint columns define separate namespaces. A metabolite called A
  and a gene called A remain two nodes and their relationship remains eligible.
  Metabolite A–gene B and metabolite B–gene A are distinct edges. Only repeated
  records of the same layer-specific pair are deduplicated.
- Bipartite MR and weight tables retain raw identifiers. ClusterONE input and
  module members use `metabolite:<URL-escaped identifier>` and
  `gene:<URL-escaped identifier>` to preserve layer identity through clustering.
  Decode the suffix with `urllib.parse.unquote` when displaying original names.

The optional legacy `all_objects` argument is supported as an endpoint coverage
check; it no longer determines node identity. The implementation ranks incident
edges directly and avoids constructing all possible node combinations.

## Calling paths and preserved settings

`corrMultiomics.py` is the only original Python caller in this checkout. Its
multiple-decay and programmatic single-decay branches now use the explicit
network type and layer-safe ClusterONE export. Both are exercised end to end.
The notebook invokes `corrMultiomics_mod.py`, which is absent in this checkout;
`setup.py` also lists that absent historical script. Those artifacts do not
establish which historical code generated any database. Use the present
`corrMultiomics.py` with the existing input arguments for future runs.

The [referenced upstream implementation](https://github.com/IntegrativeBioinformaticsLab/meantools/blob/main/mutual_ranks.py)
was inspected and has the same endpoint-concatenation/rank-overwrite algorithm.
The historical database-building version remains **unknown**. The checked-out
base commit is `4b0dcbba4a561b8d656f6649840e5a4222beca94`; provenance also
records the SHA256 of the actual modified ranking implementation.

The initial reruns reuse all stored correlations, including their signs and
existing filtering. They never read or transform an input matrix. Weighting
remains `exp(-(MR - 1) / DR)`, retaining strictly `weight > 0.01`. DR5/10/25 are
used for tomato; the available Drosophila networks use DR10/25. The bundled
ClusterONE 1.0 JAR is invoked with the original wrapper's arguments,
`java -jar <jar> <edges> -F csv`, without changing clustering parameters. Its
SHA256 is `3344c954df388f5cfe9795e0db375ae7189334a1eacb5bd3121507b33765a16b`.

Both stored and corrected weights are additionally clustered with the same
JAR, Java runtime and sorted canonical edge order. These baseline reruns help
separate rank effects from historical order/runtime differences. Historical
ClusterONE settings cannot be proven from the stored tables. The supplemental
comparison verifies that the requested cutoff reproduces the complete
historical retained edge set from its stored full MR table.

MR-producing CLI runs refuse an existing output database. Intermediate files
use a private temporary directory. Every new database includes a provenance
table. The completed audit additionally recorded source hashes, JAR/runtime
details, complete commands and source-code snapshots. A failed ClusterONE
process now raises instead of allowing a misleading successful run.

## Tomato input audit

All distinct non-self pairs occur in both directions, with **exactly equal**
correlation scores. There are no repeated directed pairs or conflicting scores.

| Network | Stored correlation rows | Self rows | Distinct non-self pairs |
| --- | ---: | ---: | ---: |
| Coexpression | 14,201,724 | 3,978 | 7,098,873 |
| Coabundance | 1,115,064 | 1,112 | 556,976 |

All 7,098,873 coexpression MRs and all 556,976 coabundance MRs decrease. Median
stored/corrected MR ratios are 2.001238 and 2.004616 respectively; ratios range
up to 4 because removing self-correlations matters most for the top ranks.

## Tomato edge retention and connectivity

Counts include every input node, including isolated nodes. Components include
singleton isolates. No previously retained edge is lost at any decay rate.

| Network | DR | Stored edges | Corrected edges | Components before → after | Isolates before → after |
| --- | ---: | ---: | ---: | ---: | ---: |
| Coexpression | 5 | 13,227 | 31,933 | 177 → 21 | 161 → 20 |
| Coexpression | 10 | 29,138 | 64,334 | 36 → 1 | 34 → 0 |
| Coexpression | 25 | 77,901 | 165,328 | 1 → 1 | 0 → 0 |
| Coabundance | 5 | 4,054 | 9,198 | 34 → 1 | 27 → 0 |
| Coabundance | 10 | 8,444 | 18,451 | 5 → 1 | 4 → 0 |
| Coabundance | 25 | 22,437 | 49,078 | 1 → 1 | 0 → 0 |

For the coexpression edges retained both before and after at DR25, the median
weight increase is 3.594485-fold. See each run's `rank_weight_comparison.json`
for whole-network MR quantiles and retained-edge weight differences.

## Tomato modules and pathway genes

ClusterONE modules overlap. Compare actual member sets rather than cluster IDs,
which are local labels that can change between runs. Exact member-set matches
and each old module's best Jaccard match are exported. Best matches are
directional, may be many-to-one, and describe membership similarity rather
than a one-to-one module assignment.

| Network | DR | Stored modules | Baseline rerun modules | Corrected modules | Stored/corrected identical member sets |
| --- | ---: | ---: | ---: | ---: | ---: |
| Coexpression | 5 | 217 | 217 | 778 | 12 |
| Coexpression | 10 | 814 | 814 | 630 | 20 |
| Coexpression | 25 | 627 | 629 | 370 | 5 |
| Coabundance | 5 | 153 | 153 | 183 | 11 |
| Coabundance | 10 | 214 | 216 | 125 | 16 |
| Coabundance | 25 | 124 | 123 | 65 | 5 |

For coexpression, baseline reruns exactly reproduce all 217 stored DR5 member
sets, 808 of 814 at DR10, and 616 of 627 at DR25. The baseline differences
preclude attributing every historical module change exclusively to the rank
fix; the matched current-runtime comparisons are also available.

All six corrected pathway-pair MRs and DR25 weights reproduce the independent
tomato audit. For example:

| Pair suffixes | Stored MR | Corrected MR | Stored DR25 weight | Corrected DR25 weight |
| --- | ---: | ---: | ---: | ---: |
| 240–260 | 12 | 4.898979 | 0.644036 | 0.855594 |
| 250–260 | 12 | 5 | 0.644036 | 0.852144 |
| 260–270 | 115.325626 | 56.284989 | 0.010327 | 0.109547 |

The full names are `Solyc12g100240.1`, `Solyc12g100250.1`,
`Solyc12g100260.1` and `Solyc12g100270.1`. No ranking or clustering code uses
these targets; the audit looked them up only after clustering.

The four genes do **not** share a corrected module at DR5, DR10 or DR25.
At DR25, corrected cluster 190 contains 240/250/260 (22 total members).
Corrected clusters 84 and 127 contain 260/270 (50 and 57 members). Thus 260
participates in overlapping modules connecting the two groups. At DR5 and
DR10, modules containing 270 remain separate from those containing 240/250/260.
These are observed clustering results, without forcing pathway membership.

## Additional available coabundance/network audit

The other coabundance database available in the adjacent MATE datasets,
`data_drosphila/drosophila_MET.sqlite`, is also affected. Its 1,129 correlation
rows contain 81 self-pairs and both orientations of 524 non-self pairs; reverse
scores agree exactly. It has 81 input nodes, including two with only diagonal
records; both remain represented as isolates in connectivity comparisons.

The adjacent Drosophila coexpression database is affected as well. Its
34,594,537 rows include 7,007 self-pairs and both orientations of 17,293,765
non-self pairs. The largest reverse-score difference is `1.4433e-15`, within
the documented tolerance. Its feature-label mapping table is copied to the new
database. Both Drosophila originals retain their source hashes.

| Network | DR | Stored edges | Corrected edges | Components before → after | Isolates before → after | Stored / baseline / corrected modules |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Drosophila coabundance | 10 | 522 | 524 | 3 → 3 | 2 → 2 | 16 / 16 / 9 |
| Drosophila coabundance | 25 | 524 | 524 | 3 → 3 | 2 → 2 | 9 / 10 / 8 |
| Drosophila coexpression | 10 | 28,562 | 79,774 | 2,144 → 295 | 2,119 → 293 | 658 / 658 / 1,174 |
| Drosophila coexpression | 25 | 99,017 | 246,850 | 463 → 1 | 459 → 0 | 1,176 / 1,176 / 1,085 |

No historical edge is lost in these comparisons. Coabundance module membership
changes even at DR25, where the retained edge set is unchanged, because the
edge weights change. Detailed member-set comparisons accompany each rerun.

## Separate normalization assessment

The original tomato RNA matrix has 3,978 features × 30 samples, and the
metabolite matrix has 1,112 features × 45 samples. Both have column means near
zero and sample standard deviations near one: they are standardized across
features within each sample. Applying additional row z-scores changes Pearson
PCCs by at most `1.56e-15`; it does not undo the sample-wise transformation.

For the unambiguous Mz reference subset (340 shared metabolite identifiers,
9 samples, 57,630 pairs), supplied columns are positive affine transformations
of author raw values. However, feature-profile correlations across samples are
different: median absolute PCC difference 0.166137, 95th percentile 0.816333,
and 13,122 pairs change sign. Membership at `abs(PCC) >= 0.1` changes for
11,415 pairs in this **subset-only diagnostic**. These are not full-network
edge or module comparisons and use neither an alternative gene matrix nor
the complete metabolite feature/sample set.

A separate normalization rerun needs the complete raw RNA inputs and an
agreed full reference/sample mapping. No normalization was changed for any
corrected network reported here. Diagnostic outputs were isolated in
`normalization_assessment/` with unchanged-input hashes and subsequently removed.

## Outputs and reproduction

The generated results under `results/mutual_rank_correction_2026-11-24/` were
removed at the user's request after validation, along with the auxiliary audit
scripts. The findings above are retained. The original output layout was:

- `coexpression/corrected.sqlite` and `coabundance/corrected.sqlite` contain
  newly computed MR, weights and ClusterONE tables using the original table
  prefixes. They also contain baseline rerun modules, comparison tables and
  provenance. Stored correlation tables are referenced by immutable source
  hashes rather than copied into these output databases.
- Each run contains `graph_comparison.csv`, `module_comparison.csv`,
  `module_best_matches.csv`, `target_pairs.csv`, `target_memberships.csv`,
  `rank_weight_comparison.json`, `provenance.json`, source snapshots, and
  ClusterONE edge files, CSV outputs and logs. Empty target CSVs are expected
  for non-tomato-gene networks.
- Additional available Drosophila networks are audited/regenerated in their
  own `drosophila_coabundance/` and `drosophila_coexpression/` directories.
- Generated SQLite databases and edge-list text files are git-ignored.
  Original external databases are opened with
  `mode=ro` and their SHA256 values are checked before and after processing.

For future network construction, use the existing correlation CLI with a
**new** output database. For example, a gene–gene network uses the same matrix
for both input arguments:

```sh
python corrMultiomics.py \
  -ft transcripts.csv -qm transcripts.csv \
  -dn coexpression_corrected.sqlite -tn coexpression \
  -mr -cl -c 0.1 -w 0.01 -mdr 5 10 25
```

For a gene–metabolite network, supply separate metabolite and transcript files
to `-ft` and `-qm`. These CLI runs construct networks from input matrices;
the auxiliary tools for the historical stored-correlation comparisons above
are no longer included.

## Validation

The full-environment suite is exercised on Python 3.8 with pandas 1.4.3.
It covers the requested numerical example, upper-triangle/symmetric inputs,
self-pair removal, repeated/shuffled edges, deterministic ties, signed negative
correlations, tolerance and conflicting duplicates, bipartite calculations
with both disjoint and overlapping identifiers, unchanged weighting/cutoff,
and unchanged inputs. The caller integration
test exercises both network types, both decay-rate branches, real ClusterONE
and refusal to overwrite an existing database, without a network-type flag.
File-identity tests cover path aliases, symlinks, hard links and distinct files
with identical contents.

The independent ranking tests are also exercised on Python 3.13 with
pandas 2.2.3 (caller integration is skipped because rdkit is not installed in
that environment). Network reruns used that environment and OpenJDK 25.0.2;
the precise versions are recorded in each provenance file.

```sh
/path/to/full-meantools-python -m unittest discover -s tests -v
```
