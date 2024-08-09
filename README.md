# TarGeneSimulationStudy

Repository containing run configurations and analysis scripts for the evaluation of semi-parametric estimators in the UK-Biobank.

## Reproducibility

Workflows were run on the Eddie High Performance Computing platform.

### Workflows

To run the workflows, the following dependencies need to be available:

- DVC >= 3.51.0
- Nextflow >= 24.04.4
- Singularity >= 3.8.6

Then run:

```bash
dvc repro
```

### Analyses

To run the analyses scripts, the following dependency needs to be installed:

- Julia >= 1.10.0

Then run:

```bash
julia --project -e'using Pkg; Pkg.instantiate()'
```

to install the project's dependencies and

```bash
julia --project analyses.jl'
```

to generate the plots.