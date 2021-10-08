## Current
* Add genbank-to-ddbj for Genbank → DDBJ conversion
* Add compare-ddbj for comparing two DDBJ annotations for evaluation
* Add a non-SO wild items to the default config for feature and qualifier names
* Fix: Respect feature-wise /transl_table value than the globally-set one.
* Fix an error when joined features have .sub_features
* Fix [source] in metadata not used when "source" feature exists in a entry

## v0.2.4 (2021-10-04)
* Expand the default config for feature and qualifier names
    * Now based on Sequence Ontology
* Create bgzip file as a new file, rather than replacing existing FASTA file

## v0.2.3 (2021-09-30)
* Switch to pysam and samtools to index FASTA (called faidx)
    * No longer creates sqlite3 database when reading FASTA
* Minor fixes

## v0.2.2 (2021-09-30)
* Rename the CLI option from --rename_setting to --config_rename
* Rename the CLI option form --filter_setting to --config_filter
* Update the feature-qualifier config for renaming
* Minor bugfixes

## v0.2.1 (2021-09-28)
* Add the CLI option --version

## v0.2.0 (2021-09-28)
* Improve memory usage when reading a large FASTA file
* Add CLI options: --rename_setting and --filter_setting to load custom settings
* Fix locations with inequalities '<' and '>' when start/stop codons are absent
* Fix a critical bug not inserting source properly
* Fix a critical bug missing subfeatures in certain cases
* Fix a critical bug ocassionally getting mRNA's ID wrong.
* Update the translate-feature renaming table
* Update the translate-feature pair filtering
* Add CLI tool: normalize-entry-names for renaming IDs in DDBJ annotation file
* Many bugfixes

## v0.1.1 (2021-09-14)
* Initial release to PyPI