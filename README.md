# GFF3-to-DDBJ
日本語版は[こちら](https://github.com/yamaton/gff3toddbj/blob/main/README-ja.md)。

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/yamaton/gff3toddbj?style=for-the-badge)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/gff3toddbj?style=for-the-badge)](https://bioconda.github.io/recipes/gff3toddbj/README.html)
[![PyPI](https://img.shields.io/pypi/v/gff3toddbj?style=for-the-badge)](https://pypi.org/project/gff3toddbj/)


## Overview

GFF3-to-DDBJ converts GFF3 and FASTA files into the [DDBJ annotation format](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#annotation) required for submission. It is the DDBJ-specific equivalent of tools like `table2asn` (NCBI) or `EMBLmyGFF3` (ENA).

View the [tests/golden](https://github.com/yamaton/gff3toddbj/tree/main/tests/golden) directory for example output (`.ann` files).



## Accuracy and Validation

Since "perfect" GFF3-to-DDBJ conversion is not formally defined, this tool uses RefSeq GFF3-GenBank correspondence as a gold standard. We validate output by:

1. Comparing `gff3-to-ddbj` results against GenBank-sourced annotations via an internal `genbank-to-ddbj` tool.
2. Passing all output through the [DDBJ BioProject/BioSample/Sequence Data (MSS) Parser](https://www.ddbj.nig.ac.jp/ddbj/parser-e.html).



## Installation

### Via Bioconda

```shell
conda create -n ddbj -c conda-forge -c bioconda gff3toddbj
conda activate ddbj
```



### Via PyPI

```
conda create -n ddbj -c conda-forge -c bioconda pip samtools
conda activate ddbj
python -m pip install gff3toddbj
```



### Via GitHub (Nightly)

```shell
conda create -n ddbj pip
conda activate ddbj
python -m pip install 'git+https://github.com/yamaton/gff3toddbj'
```



## Usage

```shell
gff3-to-ddbj \
  --fasta myfile.fa \               # Required
  --gff3 myfile.gff3 \              # Strongly Recommended (bare-minimum if absent)
  --metadata mymetadata.toml \      # Optional
  --locus_tag_prefix PREFIX_ \      # Required for BioSample
  --transl_table 1 \                # Default: 1 (Standard)
  --output output.ann               # Optional: stdout by default
```



### Argument Details

- `--locus_tag_prefix`: The prefix [assigned by BioSample](https://www.ddbj.nig.ac.jp/ddbj/file-format-e.html#locus_tag).
- `--transl_table`: Genetic code index (e.g., 11 for Bacteria). See [DDBJ Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode-e.html).



## Under the Hood

GFF3-to-DDBJ processes your data through the following pipeline:

### 1. Data Preparation

- **FASTA Compression:** If the input is standard Gzip, the tool re-compresses it using `bgzip` (e.g., creating `myfile_bgzip.fa.gz`). This enables indexing and reduces memory usage; the resulting file remains compatible with standard `gzip` tools.
- **Gap Detection:** Scans FASTA sequences for `N` runs and automatically generates `assembly_gap` features.
- **Topology Handling:** If GFF3 has `Is_circular=true`, the tool inserts a `TOPOLOGY` feature and manages [origin-spanning features](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/#origin-spanning-features).



### 2. Feature & Qualifier Mapping

- **SO-to-INSDC Translation:** Maps GFF3 "types" to DDBJ "Features" based on [Sequence Ontology](http://sequenceontology.org).
    - *Example:* `transcript` (SO:0000673) is translated to a `misc_RNA` feature.
- **Qualifier Renaming:** Converts GFF3 attributes to DDBJ-compliant qualifiers based on [renaming rules](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml).
    - *Example:* `ID=foobar` becomes `/note="ID:foobar"`.
- **Genetic Code Assignment:** Automatically adds the `/transl_table` qualifier to every `CDS` feature based on the user-provided index (default: 1).



### 3. Coordinate Processing

- **Joining:** Features sharing a parent are merged using `join()` notation. This applies to `CDS`, `exon`, `mat_peptide`, `V_segment`, `C_region`, `D-loop`, and `misc_feature`.
- **RNA/Exon Logic:** The location of joined `exons` is assigned to the parent RNA's location, and individual `exon` entries are discarded.
    - *Note:* `exons` are **not** joined if their direct parent is a `gene`.
- **Partialness:** Adds partial indicators (`<` or `>`) to `CDS` locations if start or stop codons are missing. (See: [Offset of the frame at translation initiation by codon_start](https://www.ddbj.nig.ac.jp/ddbj/cds-e.html#frame)).



### 4. DDBJ Compliance Logic (Product & Gene)

- **Product Enforcement:** To conform to [DDBJ instructions](https://www.ddbj.nig.ac.jp/ddbj/qualifiers-e.html#product), each `CDS` is restricted to a single `/product`:
    * Even if there are multiple general names for the same product, do not enter multiple names in 'product'. Do not use needless symbolic letters as delimiter for multiple names. If you would like to describe more than two names, please enter one of the most representative name in /product qualifier, and other(s) in /note qualifier.
    * If the name and function are not known, we recommend to describe as "hypothetical protein".
- **Gene Consistency:**
    - Ensures the `/gene` qualifier has a single value; additional values move to `/gene_synonym`. (Reference: [Definition of Qualifier key: /gene](https://www.ddbj.nig.ac.jp/ddbj/qualifiers-e.html#gene)).
    - Copies `/gene` and `/gene_synonym` qualifiers from parent `gene` features to all children (e.g., `mRNA`, `CDS`).



### 5. Metadata & Filtering

- **Metadata Injection:** Inserts `source` information and global qualifiers from the metadata file. See "Metadata Configuration" in "Customization" below.
- **Compliance Filtering:** Removes features and qualifiers violating the [DDBJ usage matrix](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf).
    - *Note:* The `gene` feature is discarded by default in this process.
- **Deduplication:** Removes redundant qualifier values generated during processing.



### 6. Final Formatting

- **Sorting:** Lines are ordered by start position, [feature priority](https://github.com/yamaton/gff3toddbj/blob/1cea725cca2a8f3edb45bac45d7983e255285d5e/gff3toddbj/transforms.py#L763) (placing `source` and `TOPOLOGY` at the top), and end position.

- **Validation Logs:** Displays all discarded items via `stderr`:

    ```plaintext
    WARNING: [Discarded] feature -------> gene (count: 49911)
    WARNING: [Discarded] (Feature, Qualifier) = (mRNA, Parent) (count: 57304)
    ```

    



## Customization

### Metadata Configuration

Use a TOML file (e.g., `metadata.toml`) to provide information absent from GFF3/FASTA files, such as submitter details and common qualifiers.

- **Example:** See [metadata_ddbj_example.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_ddbj_example.toml).
- **Default:** If `--metadata` is omitted, the tool uses this [default configuration](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/metadata_without_COMMON.toml).



#### Key Sections

1. **COMMON Entry**: Define `SUBMITTER`, `REFERENCE`, and `COMMENT` blocks.

2. **Global Qualifiers (DDBJ-side injection)**: Use the `[COMMON.feature]` syntax to instruct the DDBJ system to insert qualifiers into every occurrence of a feature.

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

3. **Local Injection (Tool-side injection)**: Use the `[feature]` syntax (without the `COMMON` prefix) to have `gff3-to-ddbj` explicitly insert these qualifiers into the generated `.ann` file.

    ```toml
    [assembly_gap]
    estimated_length = "<COMPUTE>"  # Automatically calculate gap size from "N" runs
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

    *Note: Currently, only `[source]` and `[assembly_gap]` are supported for local injection.*



### [Advanced] Feature and Qualifier Renaming

GFF3 and DDBJ formats do not share a 1:1 nomenclature. GFF3 "types" (column 3) map to DDBJ "Features," while GFF3 "attributes" (column 9) map to DDBJ "Qualifiers."

`gff3-to-ddbj` uses a [default translation table](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) to handle these conversions. You can override these rules using `--config_rename <FILE>`.



#### Customization Examples:

* **Renaming Types:** Map a GFF3 type to a specific DDBJ feature key.

    ```toml
    [five_prime_UTR]
    feature_key = "5'UTR"
    ```

* **Renaming Attributes:** Map GFF3 attributes to DDBJ qualifiers. Use `__ANY__` to apply a rule across all feature types.

    ```toml
    [__ANY__.ID]
    qualifier_key = "note"
    qualifier_value_prefix = "ID:"  # optional
    ```

* **Complex Translations:** Map a GFF3 type to a DDBJ feature/qualifier pair (e.g., `snRNA` to `ncRNA` with a class).

    ```toml
    [snRNA]
    feature_key = "ncRNA"
    qualifier_key = "ncRNA_class"
    qualifier_value = "snRNA"
    ```

* **Attribute-to-Feature Mapping:** Convert specific attribute values into distinct DDBJ features (e.g., `RNA` type with `biotype=misc_RNA` attribute becomes a `misc_RNA` feature).

    ```toml
    [RNA.biotype.misc_RNA]
    feature_key = "misc_RNA"
    ```



### [Advanced] Feature and Qualifier Filtering

To comply with the [DDBJ usage matrix](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf), output is filtered by a [default configuration](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/ddbj_filter.toml). Only features and qualifiers explicitly allowed in this TOML file will appear in the final output.

To use a custom filter, provide a TOML file via `--config_filter <FILE>` using the following structure:

```toml
# Only these qualifiers will be kept for the CDS feature
CDS = ["EC_number", "inference", "locus_tag", "note", "product"]
```




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



## Known Issues

### Biological & Sequence Logic

- **Trans-splicing:** The tool does not currently support coordinate correction or the `join()` syntax for features containing the `/trans_splicing` qualifier.
- **Translation Exceptions:** Coordinate handling for `/transl_except` at start or stop codons is not yet implemented.
- **Missing Qualifiers:** The tool does not automatically generate a `/translation` qualifier when an `/exception` qualifier is present, which may lead to DDBJ validation errors.
- **Inter-base Coordinates:** "Between-position" locations (e.g., `123^124`) are not currently supported and may be processed incorrectly.

### Performance

- **Execution Speed:** To ensure maximum accuracy, the tool currently utilizes a single-process architecture. Expect longer runtimes on large genomic datasets.



## Acknowledgments

The design of GFF3-to-DDBJ is inspired by [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3), a versatile tool used for converting GFF3 data into the EMBL annotation format.
