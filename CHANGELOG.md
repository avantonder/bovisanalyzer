# avantonder/bovisanalyzer: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.1 - [22/11/22]

- Fix Kraken 2 database bugs

## v1.2 - [14/11/22]

- Add test profile and dataset so pipeline can be run with test.conf
- Rewrite python scripts
- Add skip_alignment, skip_kraken2 and skip_clusters parameters
- Update nextflow_schema.json
- Add parameters.md to docs

## v1.1 - [26/10/22]

- Fix masking bug
- Output Picard MARKDUPLICATES metrics for multiQC
- Delete unnecessary subworkflows and modules
- Correct samtools index version from 1.14 to 1.15.1
- Edit documentation
- Standardise old modules to current syntax
- Add skip_alignment option to make alignment creation optional
- Edit modules.config to output low_quality_pseudogenomes.tsv
- Update spoligotype database

## v1.0 - [21/10/22]

First release version of avantonder/bovisanalyzer.

## v1.0dev - [28/04/22]

Initial release of avantonder/bovisanalyzer, created with the [nf-core](https://nf-co.re/) template.

