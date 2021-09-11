#!/usr/bin/env python

import subprocess
import pathlib
import uuid
import os

WORK_DIR = "tests/golden/"

class BunchOfFiles(object):
    def __init__(self, reference, input_gff3, input_fasta, input_metadata, transl_table=1):
        self.ref = reference
        self.gff3 = input_gff3
        self.fasta = input_fasta
        self.metadata = input_metadata
        self.transl_table = transl_table


def get_command(data: BunchOfFiles, output: str):

    command = """\
    python src/main.py \
        --gff3={gff3} \
        --fasta={fasta} \
        --transl_table={transl_table} \
        --metadata={metadata} \
        --output={output} \
        2> /dev/null
    """.format(
        gff3=data.gff3,
        fasta=data.fasta,
        metadata=data.metadata,
        transl_table=data.transl_table,
        output=output,
    )

    return command


def _runner(data: BunchOfFiles):
    output_file = pathlib.Path(WORK_DIR) / (str(uuid.uuid1()) + ".ann")
    cmd = get_command(data, output_file)
    subprocess.run(cmd, shell=True)

    with open(data.ref, "r") as f_ref, open(output_file, "r") as f_out:
        content_ref = f_ref.read()
        content_out = f_out.read()

    assert content_ref == content_out

    if output_file.exists():
        os.remove(output_file)
    else:
        raise FileNotFoundError("Failed to create output file in golden test!; {}".format(output_file.as_posix))

def test_golden():
    testdata0 = BunchOfFiles(
        WORK_DIR + "output_GCF_000280675.1_ASM28067v1_genomic.ann",
        WORK_DIR + "GCF_000280675.1_ASM28067v1_genomic.gff.gz",
        WORK_DIR + "GCF_000280675.1_ASM28067v1_genomic.fna.gz",
        WORK_DIR + "metadata.toml",
        4,  # transl_table
    )
    _runner(testdata0)
