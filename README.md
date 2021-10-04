# GFF3-to-DDBJ
日本語版は[こちら](https://github.com/yamaton/gff3toddbj/blob/main/README-ja.md)。

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/yamaton/gff3toddbj?style=for-the-badge)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/gff3toddbj?style=for-the-badge)](https://bioconda.github.io/recipes/gff3toddbj/README.html)
[![PyPI](https://img.shields.io/pypi/v/gff3toddbj?style=for-the-badge)](https://pypi.org/project/gff3toddbj/)


[TOC]

## Table of Contents
* [What is this?](#what-is-this)
* [Initial setup](#initial-setup)
  + [Install with bioconda](#install-with-bioconda)
  + [Install with pip](#install-with-pip)
  + [Install from the source](#install-from-the-source)
* [Create DDBJ annotation from GFF3 and FASTA](#create-ddbj-annotation-from-gff3-and-fasta)
  + [Run `gff3-to-ddbj`](#run-gff3-to-ddbj)
* [Under the Hood](#under-the-hood)
* [Customize the behavior](#customize-the-behavior)
  + [Metadata file](#metadata-file)
  + [[Advanced] Rename features and qualifiers](#advanced-rename-features-and-qualifiers)
    + [Rename types/feature keys](#rename-typesfeature-keys)
    + [Rename attributes/qualifier keys](#rename-attributesqualifier-keys)
    + [Translate GFF3 types to features with qualifiers](#translate-gff3-types-to-features-with-qualifiers)
    + [Translate (type, attribute) items to features](#translate-type-attribute-item-to-feature)
    + [Run with custom configuration](#run-with-custom-configuration)
  + [[Advanced] Filter features and qualifiers](#advanced-filter-features-and-qualifiers)
* [Troubleshooting](#troubleshooting)
  + [Validate GFF3](#validate-gff3)
  + [Split FASTA from GFF3 (if needed)](#split-fasta-from-gff3-if-needed)
  + [Normalize entry names (if needed)](#normalize-entry-names-if-needed)
* [Credit](#credit)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## What is this?

GFF3-to-DDBJ creates [the annotation file for submission to DDBJ](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#annotation) by taking GFF3 and FASTA files as input. It also works with FASTA alone.

Analogous programs are [GAG](https://github.com/genomeannotation/GAG) for submissions to NCBI, and [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) for submissions to EMBL.

Please take a look at our [test directory](https://github.com/yamaton/gff3toddbj/tree/main/tests/golden) for examples. Files ending with .ann are the DDBJ annotation files produced by thie program.


## Initial setup

### Install with bioconda

```shell
# Create a conda environment named "ddbj", and install relevant packages from bioconda channel
conda create -n ddbj -c bioconda -c conda-forge gff3toddbj

# Activate the environment "ddbj"
conda activate ddbj
```

### Install with pip

```shell
# Create a conda environment named "ddbj" and install pip
conda create -n ddbj pip

# Activate the environment "ddbj"
conda activate ddbj

# Need bgzip executable in samtools
conda install -c bioconda samtools

# Install from pip
pip install gff3toddbj
```


### Install from the source

```shell
# Download
wget https://github.com/yamaton/gff3_to_ddbj/archive/refs/heads/main.zip

# Extract, rename, and change directory
unzip main.zip && mv gff3toddbj-main gff3toddbj && cd gff3toddbj

# Create a conda environment named "ddbj"
conda create -n ddbj

# Activate the environment "ddbj"
conda activate ddbj

# Install dependencies to "ddbj"
conda install -c bioconda -c conda-forge biopython bcbio-gff toml setuptools pysam samtools

# Install gff3-to-ddbj and extra tools
python setup.py install
```



## Create DDBJ annotation from GFF3 and FASTA



### Run `gff3-to-ddbj`

Let's run the main program to get some ideas.

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \                # bare-minimum output if omitted
  --fasta myfile.fa \                 # <<REQUIRED>>
  --metadata mymetadata.toml \        # example metadata used if omitted
  --locus_tag_prefix MYOWNPREFIX_ \   # default is "LOCUSTAGPREFIX_"
  --transl_table 1 \                  # default is 1
  --output myawesome_output.ann       # standard output if omitted
```

Here is the options:
* `--gff3 <FILE>` takes GFF3 file
* `--fasta <FILE>` takes FASTA file
* `--metadata <FILE>` takes the metadata file in TOML
* `--locus_tag_prefix <STRING>` takes the prefix of locus tag [obtained from BioSample](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#locus_tag). You can skip this for now.
* `--transl_table <INT>`: Choose appropriate one from [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode-e.html). The default value is 1 ("standard").
* `--output <FILE>` sets the path the annotation output.



## Under the Hood

Here is the list of operations `gff3-to-ddbj` will do:

* Re-compress FASTA with [bgzip](https://www.htslib.org/doc/bgzip.html) if the FASTA input is compressed with gzip
  * A bgzip file is created if absent like `myfile_bgzip.fa.gz`.
  * For indexing and saving memory
  * The bgzip file should be compatible with gzip

* Rename features and qualifiers following the [renaming scheme](#advanced-rename-features-and-qualifiers).

* Search for `assembly_gap`s in FASTA.

* Add `/transl_table` to each CDS.

* Insert `source` information from the [metadata fie](#metadata-file).

* Merge `CDS`s having the same parent with `join` notation.

* Merge `mRNA` and `exon` in GFF3 and create `mRNA` feature with `join` notation.

* Modify locations with inequality signs (`<` and `>`) if start/stop codon is absent.
  * See [Offset of the frame at translation initiation by codon_start](https://www.ddbj.nig.ac.jp/ddbj/cds-e.html#frame)

* Let CDS have a single `/product` value: Set it to "hypothetical protein" if absent. Move the rest of exising values to `/note`.

  * This is to conform the [instruction](https://www.ddbj.nig.ac.jp/ddbj/qualifiers-e.html#product) on `/product`.

    > * Even if there are multiple general names for the same product, do  not enter multiple names in 'product'. Do not use needless symbolic  letters as delimiter for multiple names. If you would like to describe  more than two names, please enter one of the most representative name in /product qualifier, and other(s) in /[note](https://www.ddbj.nig.ac.jp/ddbj/qualifiers-e.html#note) qualifier.
    >
    > * If the name and function are not known, we recommend to describe as "hypothetical protein".

* Remove duplicates in qualifier values.

* Sort lines in annotation

* Filter features and qualifiers following [the matrix](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf).



## Customize the behavior

### Metadata file

To enter information missing in GFF3 or FASTA, such as submitter names and certain qualifier values, you need to feed a metadata file in TOML, say `mymetadata.toml`. Take a look at [an example](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_ddbj_example.toml) matching [the example annotation in the DDBJ page](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#annotation).

The file accommodates following and they are all optional. That is, GFF3-to-DDBJ works even with an empty file.

* Basic features in the [COMMON](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#common) entry

  * ... such as `SUBMITTER`, `REFERENCE`, and `COMMENT`.

* "meta-description" in the COMMON entry

  * Here is an example with this notation:

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

  * DDBJ annotation supports "meta" values with features under COMMON such that the items are inserted to each occurrence **in the resulting flat file** produced by DDBJ. Here is an example to insert `assembly_gap` feature under `COMMON` entry.

* Feature-qualifier items inserted to each occurrence

  * Here is an example: Difference from the previous case is only at `[assembly_gap]` as opposed to`[COMMON.assembly_gap]`.

    ```toml
    [assembly_gap]
    estimated_length = "unknown"   # Set it "<COMPUTE>" to count the number of N's
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

  * While this should work effectively the same as the "meta-description" item above, use this notation if you insert values repeatedly **in the annotation file** produced by GFF3-to-DDBJ.

  * Currently supporting `[source]` and `[assembly_gap]` only.

  * If both `[COMMON.assembly_gap]` and `[COMMON.assembly]` exist in the metadata file, `gff3-to-ddbj` takes the one with COMMON.

For more examples, see [WGS in COMMON](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=1110334278) and [WGS](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=382116224) provided by DDBJ as annotation examples, and corresponding metadata files [metadata_WGS_COMMON.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS_COMMON.toml) and [metadata_WGS.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS.toml) in this repository.


### [Advanced] Rename Features and Qualifiers

GFF3 and DDBJ annotation have rough correspondence like:

1. GFF3 column 3 "type" →  DDBJ annotation column 2 as "Feature"
2. GFF3 column 9 "attribute" →  DDBJ annotation column 4 and 5 as "Qualifier key", and "Qualifier value"

but nomenclatures in GFF3 often do not conform the annotations set by INSDC. Furthermore, DDBJ lists up the [feature-qualifier pairs they accepts](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit#gid=0), a subset of the INSDC definitions.

To meet convensions with the requirement, GFF3-to-DDBJ comes with [a default configuration in TOML](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) to rename (or even translate) feature keys and qualifier keys/values. Note that [the Sequence Ontology](http://sequenceontology.org/browser) is helpful in translating a type into a INSDC feature and qualifier(s).

Here is how to customize the renaming configuration.

#### Rename types/feature keys

The default setting renames `five_prime_UTR` "type" in GFF3 into `5'UTR` "feature key" in the annotation. This transformation is expressed in TOML as follows:

```toml
[five_prime_UTR]
feature_key = "5'UTR"
```

#### Rename attributes/qualifier keys

This is about renaming attributes under arbitrary types. By default, `ID=foobar` "attribute" in a GFF3 becomes `/note="ID:foobar"` qualifier in the annotation. (Here I follow the convention putting slash (like `/note`) to denote qualifier. But DDBJ annotation does NOT include slash hence no slash is used in any of TOML files.)

Here is the TOML defining the transformation. `__ANY__` is the special name representing arbitrary types. `ID` is the original attribute key. `note` is the name of corresponding qualifier key. `ID:` is attached as the prefix of the qualifier value.

```toml
[__ANY__.ID]
qualifier_key = "note"
qualifier_value_prefix = "ID:"  # optional
```

#### Translate GFF3 types to features with qualifiers

Sometimes we want to replace a certain types with features WITH qualifiers. For example, `snRNA` is an invalid feature in INSDC/DDBJ hence we replace it with `ncRNA` feature with `/ncRNA_class="snRNA"` qualifier. Such transformation is written in TOML as following.

```toml
[snRNA]
feature_key = "ncRNA"
qualifier_key = "ncRNA_class"
qualifier_value = "snRNA"
```

#### Translate (type, attribute) items to features

Example: some annotation programs produce a GFF3 line containing `RNA` as the type and `biotype=misc_RNA` as one of the attributes. Then it should be translated to `misc_RNA` feature in annoation.

```toml
[RNA.biotype]
attribute_value = "misc_RNA"
feature_key = "misc_RNA"
```


#### Run with custom configuration

To feed a custom translation table, use the CLI option:

* `--config_rename <FILE>`

And here is an example call:

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata mymetadata.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --config_rename my_translate_features_qualifiers.toml \  # Set your customized file here
  --output myawesome_output.ann
```

### [Advanced] Filter features and qualifiers

DDBJ specifies recommended [Feature/Qualifier usage matrix](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf). To conform this rule, features and qualifiers appearing in the annotation output are filtered by [the filtering file in TOML](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/ddbj_filter.toml) by default. The file is in TOML format with the structure like this:

```toml
CDS = [
"EC_number",
"inference",
"locus_tag",
"note",
"product",
]

exon = [
"gene",
"locus_tag",
"note",
]
```

The left-hand side of the equal sign `=` represents an allowed feature key, and the right-hand side is a list of allowed qualifier keys. In this example, only `CDS` and `exon` features will show up in the annotation, and qualifiers are limited to the listed items. To customize this filtering function, edit the TOML file first and pass the file with the CLI option:

* `--config_filter <FILE>`


## Troubleshooting

### Validate GFF3

It might be a good practice to validate your GFF3 files. [GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) is useful though the file size is limited to 50MB.



### Split FASTA from GFF3 (if needed)

GFF3_to_DDBJ does not work when GFF3 contains FASTA information inside with `##FASTA` directive. Attached tool `split-fasta` reads a GFF3 file and saves GFF3 (without FASTA info) and FASTA.

```shell
split-fasta path/to/myfile.gff3 --suffix "_splitted"
```

This creates two files, `myfile_splitted.gff3` and `myfile_splitted.fa`.



### Normalize entry names (if needed)

Letters like `=|>" []` are not allowed in the 1st column (= "Entry") of the DDBJ annotation.  The attached program `normalize-entry-names` renames such entries. This program converts an ID like `ERS324955|SC|contig000013` into `ERS324955:SC:contig000013` for example.

```shell
normalize-entry-names myannotation_output.txt
```

This command create as files `myannotation_output_renamed.txt` *if* the invalid letters are found. Otherwise, you'll see no output.





## Credit

GFF3-to-DDBJ's design is deeply indebted to [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3), a versatile coversion for EMBL annotation format.
