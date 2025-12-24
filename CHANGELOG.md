## v0.4.3 (2025-12-24)
* fix: Create SeqRecord with Seq, not with a plain string for Biopython 1.86

## v0.4.2 (2025-12-24)
* fix: restore compatibility with Biopython 1.80+ unhashable SeqFeatures
* Pin biopython 1.81 for compatibility

## v0.4.1 (2025-12-23)
* Pin dependencies to make existing code work.
* Update to modern Python packaging (PEP 621)

## v0.4.0 (2021-10-22)
* Change the structure of renaming config in TOML.
    * We no longer use "attribute_value = ..."
* Simplify internals of `genbank-to-ddbj`.
* Change `compare-ddbj` to always count items as multiset elements.
* Accept case variants in feature/qualifier names.
* Add TOPOLOGY feature if `Is_circular=true` exists in `source` feature
    * Add handling of origin-spanning features in circular genome
* Merge exons with their parent if its type is `*_segment`, `*_region`, or `*RNA*`.
* Add new renaming rule: Replace a certain qualifier (key, value) with another
    * This enables the replacement: `/exception="ribosomal slippage"` --> `/ribosomal_slippage`.
* Add the renaming /genome="mitochondrion" -> /organelle="mitochondrion" to the default config
* Enforce single value to `/gene` if it has multiple. Put the rest to `/gene_synonym`.
* Do not add partial markups (`<` or `>`) if `/pseudo` or `/pseudogene` exist.
* Change start codons to ATG only when the genetic code is 1 (`/transl_table=1`).
* Do not add `assembly_gap` if the gap is < 10bp long.
* Join features if `/ribosomal_slippage` exists.
* Copy qualifiers `/gene` and `/gene_synonym` to children if a `gene` feature has them.
* Add an renaming to `/organelle`, `/host`, `/sub_strain`, `/plasmid` to the default config
* Fix a bug in DDBJ parser used in evaluation
    * It gave wrong position when location string is digits.
* Add preliminary handling of "between-position" location like `138683^138684`
    * Without `BetweenPosition` because `Bio.SeqIO.parse()` creates `FeatureLocation` instead.
* Fix representation of a single-position location with a partial markup (form `<100..100,200..300` to `<100,200..300`).
* Fix missing FASTA entries that are absent in GFF3
    * Respect the entry order in GFF3 (if provided), and append the only-in-FASTS entries.


## v0.3.0 (2021-10-12)
* Add `genbank-to-ddbj` for Genbank → DDBJ conversion
* Add `compare-ddbj`s for comparing two DDBJ annotations for evaluation
* Add a non-SO wild items to the default config for feature and qualifier names
* Update the default config: replace hyphens in qualifier keys with underscores.
* Fix: Respect feature-wise /transl_table value than the globally-set one.
* Fix an error when joined features have .sub_features
* Fix [source] in metadata not used when "source" feature exists in a entry
* Change not to join features directly under "gene"
* Fix the default config for feature and qualifier names around /pseudo
* Fix a bug getting incorrect codon when CompoundLocation has a part < 3bp.
* Fix start codon detection: no longer change /codon_start value even when !=1
* Fix stop codon detection
    * consider stop codon does not exists whenever the length
      (after subtraction by "phase") is not a multiple of 3.

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