# Snakemake workflow: snakemake_learning_rnaseq

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/snakemake_learning_rnaseq.svg?branch=master)](https://travis-ci.org/snakemake-workflows/snakemake_learning_rnaseq)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

* chaodi (@dic)

## the snakemake template from the tutorial
```
├── config
│   ├── config.yaml
│   └── samples.tsv
├── LICENSE
├── README.md
├── resources
│   └── README.md
├── results
│   ├── logs
│   ├── plots
│   ├── README.md
│   └── tables
└── workflow
    ├── envs
    │   └── myenv.yaml
    ├── report
    │   ├── some-plot.rst
    │   └── workflow.rst
    ├── rules
    │   ├── common.smk
    │   └── other.smk
    ├── schemas
    │   ├── config.schema.yaml
    │   └── samples.schema.yaml
    ├── scripts
    │   └── common.py
    └── Snakefile
```
