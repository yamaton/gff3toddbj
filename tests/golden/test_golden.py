#!/usr/bin/env python

from typing import Union
import subprocess
import pathlib
import uuid

WORK_DIR = pathlib.Path(__file__).parent
TEST_OUTPUT_DEFAULT_PREFIX = "testout_"


# [TODO] Rewrite with dataclass when support for Python 3.6 is gone.
class BunchOfFiles(object):
    def __init__(
        self, reference, input_gff3, input_fasta, input_metadata, transl_table=1
    ):
        self.ref = reference
        self.gff3 = input_gff3
        self.fasta = input_fasta
        self.metadata = input_metadata
        self.transl_table = transl_table


def get_command(data: BunchOfFiles, output: Union[str, pathlib.Path]):
    command = """\
    gff3-to-ddbj \
        --gff3={gff3} \
        --fasta={fasta} \
        --transl_table={transl_table} \
        --metadata={metadata} \
        --output={output} \
        2> stderr.log
    """.format(
        gff3=data.gff3,
        fasta=data.fasta,
        metadata=data.metadata,
        transl_table=data.transl_table,
        output=output,
    )
    return command


def runner(data: BunchOfFiles, prefix:str=""):
    filename = "_".join([TEST_OUTPUT_DEFAULT_PREFIX, prefix, str(uuid.uuid1()) + ".ann"])
    output_file = WORK_DIR / filename
    cmd = get_command(data, output_file)
    subprocess.run(cmd, shell=True)

    with open(data.ref, "r") as f_ref, open(output_file, "r") as f_out:
        content_ref = f_ref.read()
        content_out = f_out.read()

    assert content_ref == content_out
    output_file.unlink(missing_ok=True)


def tests_refseq():
    ## Testing against GCF_000280675.1_ASM28067v1
    testdata_refseq = BunchOfFiles(
        WORK_DIR / "expected_GCF_000280675.1_ASM28067v1_genomic.ann",
        WORK_DIR / "GCF_000280675.1_ASM28067v1_genomic.gff.gz",
        WORK_DIR / "GCF_000280675.1_ASM28067v1_genomic.fna.gz",
        WORK_DIR / "metadata.toml",
        4,  # transl_table
    )
    runner(testdata_refseq, "refseq")


def tests_augustus():
    testdata_augustus = BunchOfFiles(
        WORK_DIR / "augustus.ann",
        WORK_DIR / "augustus.gff3",
        WORK_DIR / "augustus.fa",
        WORK_DIR / "metadata_without_COMMON.toml",
        1,  # transl_table
    )
    runner(testdata_augustus, "augustus")


def tests_maker():
    testdata_maker = BunchOfFiles(
        WORK_DIR / "maker.ann",
        WORK_DIR / "maker.gff3",
        WORK_DIR / "maker.fa",
        WORK_DIR / "metadata_without_COMMON.toml",
        1,  # transl_table
    )
    runner(testdata_maker, "maker")


def tests_prokka():
    testdata_prokka = BunchOfFiles(
        WORK_DIR / "prokka.ann",
        WORK_DIR / "prokka.gff3",
        WORK_DIR / "prokka.fa",
        WORK_DIR / "metadata_without_COMMON.toml",
        11,  # transl_table
    )
    runner(testdata_prokka, "prokka")
