General documentation for MEANtools
===================================

Read the documentation here:

https://meantools.readthedocs.io/en/latest/

Maintainer
==========
Dr. Kumar Saurabh Singh (kumarsaurabh.singh@maastrichtuniversity.nl)


Principal Investigator
======================
Dr. Marnix Medema (marnix.medema@wur.nl)




MEANtools is licensed under the [MIT License](LICENSE.md)

# Mutual-rank construction

MR runs automatically use a unimodal network when `-ft` and `-qm` refer to the
same input file, and a bipartite network when they refer to different files.
For gene–gene or metabolite–metabolite networks, pass the same matrix path to
both existing input arguments. Separate files keep separate node layers even
when their feature names or contents match. No network-type option is needed.
Use a new output database. Bipartite module members carry layer prefixes; MR
and weight tables retain the original identifiers.

See [ranking conventions, rerun commands, and validation](docs/mutual_rank_correction.md).

# HPC execution

The `HPC/` directory provides a single-node SLURM implementation for larger
correlation and mutual-rank analyses. It follows the same filtering, automatic
network selection, corrected ranking, weighting and ClusterONE conventions as
the standard `corrMultiomics.py` runner.

| Script | Purpose |
| --- | --- |
| [corrMultiomics_HPC.py](HPC/corrMultiomics_HPC.py) | Computes correlations in blocks using a bounded thread pool, streams results into SQLite, and calculates mutual ranks using scratch-disk or RAM storage. |
| [corrMultiomics.sh](HPC/corrMultiomics.sh) | Submits the Python runner through SLURM, configures CPU/thread limits, creates job-local scratch space, and gives each job a separate output database. |

Keep `HPC/` inside the MEANtools checkout so the runner can import the shared
`mutual_ranks.py`. The Python dependencies are NumPy, SciPy, pandas and
`threadpoolctl`; disk ranking requires SQLite >=3.25. ClusterONE requires Java
and the bundled `cluster_one-1.0.jar`.

Before submitting, edit the paths, environment activation, account, partition,
CPU, memory and time settings in `HPC/corrMultiomics.sh` for your cluster:

```bash
# Default: store mutual-rank intermediates on scratch disk.
sbatch HPC/corrMultiomics.sh

# Alternatively, use RAM with an estimated working-memory budget of 20 GiB.
sbatch --export=ALL,MR_STORAGE=ram,MR_MEMORY_GB=20 HPC/corrMultiomics.sh
```

The RAM budget must fit within the job's memory allocation; it does not request
additional memory from SLURM. Both modes write results to SQLite and require a
new output database for MR construction. As with the standard runner, use the
same matrix file for both `-ft` and `-qm` for a unimodal network, or separate
input files for a bipartite network.

See the [HPC documentation](HPC/README.md) for direct Python usage, resource
controls, output tables and validation details.
