# gff3_to_ddbj
Convert GFF3 + FASTA to DDZJ annotation

## How to use

1. Edit metadata in `metadata.toml` for the `COMMON` entry, `source` features, and `assembly_gap` features.
2. Run

```shell
python main.py --gff3 <path-to-GFF3> --fasta <path-to-FASTA> --locus_tag_prefix MYLOCUSTAG_ > ddbj_annotation.txt
```
