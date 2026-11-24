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

# HPC implementation

The `HPC/` scripts add block-based correlation calculations, controlled
threading, streamed SQLite output, and disk or RAM storage for mutual ranks.
These changes reduce per-pair processing overhead and control intermediate
memory use, helping larger analyses run on a single HPC node. Both storage
modes follow the corrected mutual-rank conventions used by the standard runner.

The included SLURM submission script manages resource settings, temporary
scratch files and separate output databases for each job.

See the [HPC README](HPC/README.md) for implementation details, requirements
and usage examples.
