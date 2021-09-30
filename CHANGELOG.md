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