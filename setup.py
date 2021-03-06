import setuptools
import pathlib

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# sets __version__
with open(pathlib.Path("gff3toddbj") / "version.py") as f:
    exec(f.read())

setuptools.setup(
    name="gff3toddbj",
    version=__version__,
    author="Yamato Matsuoka",
    author_email="yamaton@gmail.com",
    description="Create a DDBJ annotation file from GFF3 and FASTA files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yamaton/gff3toddbj",
    project_urls={
        "Bug Tracker": "https://github.com/yamaton/gff3toddbj/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],

    packages=setuptools.find_packages(),
    package_data={"": ["*.toml"]},
    python_requires=">=3.6",
    install_requires=[
        "biopython >= 1.75",
        "bcbio-gff >= 0.6.6",
        "pysam",
        "toml",
    ],
    entry_points={"console_scripts":
        [
            "gff3-to-ddbj = gff3toddbj:main",
            "genbank-to-ddbj = gff3toddbj.genbank_to_ddbj:main",
            "split-fasta = tools.splitfasta:main",
            "normalize-entry-names = tools.normalize_seqids:main",
            "list-products = tools.list_products:main",
            "compare-ddbj = gff3toddbj.evaluate:main"
        ]},
)
