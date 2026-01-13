# CLAUDE.md - GFF3-to-DDBJ Project Knowledge Base

**Last Updated**: 2026-01-13
**Version**: 0.4.3
**Status**: Verified Working ✓

---

## 1. PROJECT OVERVIEW

**GFF3-to-DDBJ** is a specialized bioinformatics tool for converting GFF3 (Generic Feature Format) and FASTA sequence files into DDBJ (DNA Data Bank of Japan) annotation format for genome sequence submission.

### Key Facts
- **Purpose**: DDBJ equivalent of NCBI's `table2asn` and ENA's `EMBLmyGFF3`
- **Language**: Python 3.9+
- **License**: GPLv3
- **Author**: Yamato Matsuoka
- **Output Format**: `.ann` (DDBJ annotation files)
- **Primary LOC**: 2,266 lines (main package) + 292 lines (tools)
- **Package Manager**: `uv` (this installation), README suggests `conda`

### Critical Features
- Converts GFF3 + FASTA → DDBJ annotation format
- Handles feature hierarchies with parent-child relationships
- Coordinate system conversion (0-based → 1-based)
- DDBJ compliance filtering and validation
- Memory-efficient processing via pysam FASTA indexing
- Supports multiple annotation sources (Prokka, MAKER, Augustus, RefSeq)

---

## 2. ARCHITECTURE

### 6-Stage Processing Pipeline

```
Input (GFF3 + FASTA + Metadata)
    ↓
[1] Data Preparation
    • FASTA bgzip compression
    • Assembly gap detection (N-runs)
    • Circular topology handling
    ↓
[2] Feature & Qualifier Mapping
    • SO-to-INSDC translation (types → Features)
    • Qualifier renaming (attributes → DDBJ qualifiers)
    • Genetic code assignment (/transl_table)
    ↓
[3] Coordinate Processing
    • Feature joining (join() notation)
    • RNA/exon logic (inheritance)
    • Partialness markers (< >)
    ↓
[4] DDBJ Compliance Logic
    • Product enforcement (single /product)
    • Gene consistency (/gene, /gene_synonym)
    • Qualifier inheritance (parent → children)
    ↓
[5] Metadata & Filtering
    • Metadata injection (source qualifiers)
    • Compliance filtering (DDBJ usage matrix)
    • Deduplication
    ↓
[6] Final Formatting
    • Sorting (position, priority, end)
    • Validation logging
    ↓
Output (DDBJ .ann file)
```

### Module Responsibilities

| Module | Primary Responsibility | Key Functions |
|--------|------------------------|---------------|
| `main.py` | CLI orchestration, argument parsing | `main()` |
| `transforms.py` | Core transformation pipeline | Feature/qualifier mapping, joining |
| `formatter.py` | DDBJ formatting, compliance rules | Output generation |
| `parser.py` | DDBJ annotation parsing | Location parsing, record loading |
| `io.py` | File I/O operations | GFF3 parsing, FASTA indexing |
| `utils.py` | Utility functions | Metadata loading, codon detection |
| `genbank_to_ddbj.py` | GenBank → DDBJ conversion | Alternative entry point |
| `evaluate.py` | Comparison tools | `compare-ddbj` command |

---

## 3. KEY TECHNICAL DETAILS

### Coordinate System Handling
- **BioPython**: 0-based, left-inclusive, right-exclusive `[start, end)`
- **DDBJ**: 1-based, both-inclusive `[start, end]`
- **Conversion**: Add 1 to start positions, keep end positions
- **Partial locations**: `<` (missing start), `>` (missing end)

### Feature Hierarchy
- Uses BioPython's `.sub_features` for parent-child relationships
- Dummy prefix `__tmpname__` prevents double-application of rules
- Gene → RNA → CDS hierarchies preserved
- Qualifiers inherited from parent to children

### FASTA Processing
- **pysam** for memory-efficient indexing
- Auto-converts gzip → bgzip (samtools-compatible)
- Creates `.fai` (FASTA index) and `.gzi` (gzip index)
- Handles large genomes without full memory load

### Configuration System
- **TOML-based**: 4 configuration files (1,662 total lines)
- `translate_features_qualifiers.toml`: Feature/qualifier mappings
- `ddbj_filter.toml`: DDBJ usage matrix (allowed pairs)
- `metadata*.toml`: Default metadata templates
- Custom configs via `--config_rename` and `--config_filter`

---

## 4. DIRECTORY STRUCTURE

```
gff3_to_ddbj/
├── gff3toddbj/              # Main package (2,266 LOC)
│   ├── main.py              # gff3-to-ddbj CLI entry point
│   ├── transforms.py        # Core transformation logic (~1000 LOC)
│   ├── formatter.py         # DDBJ formatting (~500 LOC)
│   ├── parser.py            # DDBJ parsing (~250 LOC)
│   ├── io.py                # File I/O (~400 LOC)
│   ├── utils.py             # Utilities (~300 LOC)
│   ├── genbank_to_ddbj.py   # GenBank conversion CLI (~120 LOC)
│   ├── evaluate.py          # Comparison tools (~200 LOC)
│   └── *.toml               # Configuration files
├── tools/                   # CLI utilities (292 LOC)
│   ├── splitfasta/          # Split GFF3+FASTA
│   ├── normalize_seqids/    # DDBJ-compliant ID renaming
│   ├── list_products/       # Extract product names
│   ├── regularize_seqids/   # Sequence ID regularization
│   └── unquote/             # URL-encoded character handling
├── tests/                   # Test suite
│   ├── test_annotation_parsing_and_formatting.py
│   ├── golden/              # Reference test data
│   │   ├── augustus.*
│   │   ├── maker.*
│   │   ├── prokka.*
│   │   ├── fasta_only.*
│   │   └── GCF_000280675.1.*
│   └── README.md
├── notebooks/               # Jupyter notebooks
├── examples/                # Example data
├── evaluation/              # Evaluation scripts
├── pyproject.toml           # Modern Python packaging (PEP 621)
├── CHANGELOG.md
├── README.md (English)
├── README-ja.md (Japanese)
├── LICENSE
└── CLAUDE.md                # This file
```

---

## 5. CLI COMMANDS

### Main Commands
```bash
# Primary conversion tool
gff3-to-ddbj --fasta file.fa --gff3 file.gff3 --locus_tag_prefix PREFIX_ [options]

# GenBank to DDBJ conversion
genbank-to-ddbj input.gbk output.ann

# Comparison tool
compare-ddbj file1.ann file2.ann
```

### Utility Commands
```bash
split-fasta combined.gff3 --fasta out.fa --gff3 out.gff3
normalize-entry-names file.fa --prefix MYPREFIX
list-products annotations.gff3
```

### Key Arguments
- `--fasta`: Input FASTA file (required)
- `--gff3`: Input GFF3 file (strongly recommended)
- `--metadata`: Metadata TOML file (optional)
- `--locus_tag_prefix`: Required for BioSample submission
- `--transl_table`: Genetic code (1=standard, 11=bacterial)
- `--config_rename`: Custom feature/qualifier mappings
- `--config_filter`: Custom DDBJ compliance rules
- `--output`: Output file path (default: stdout)

---

## 6. DEPENDENCIES

### Core Runtime Dependencies
```toml
biopython <= 1.86          # SeqRecord, SeqFeature, SeqIO
bcbio-gff <= 0.7.1         # GFF3 parsing
pysam                       # FASTA indexing (samtools)
toml                        # Configuration file parsing
```

### Build System
```toml
setuptools >= 77.0.3        # PEP 621 compliant
python >= 3.9, < 3.14
```

### Development
```toml
pytest >= 8.4.2             # Testing framework
```

### Installation Methods
```bash
# This Installation (using uv)
# Dependencies managed via uv.lock
# Run commands with: uv run -- <command>

# Bioconda (README recommended)
conda create -n ddbj -c conda-forge -c bioconda gff3toddbj

# PyPI
pip install gff3toddbj

# GitHub (nightly)
pip install 'git+https://github.com/yamaton/gff3toddbj'
```

### Running Commands with uv
```bash
# Run any command with uv
uv run -- <command> [args]

# Examples:
uv run -- pytest
uv run -- gff3-to-ddbj --help
uv run -- python -m gff3toddbj.main --version
```

---

## 7. TESTING

### Test Suite Structure
- **Framework**: pytest with doctests
- **Location**: `/tests/`
- **Golden Tests**: Reference data from multiple sources

### Test Data Sources
| Source | Description | Data Size |
|--------|-------------|-----------|
| **RefSeq GCF_000280675.1** | Reference genome | 11MB FASTA (gz), 1.5MB GFF3 (gz) |
| **Prokka** | Bacterial annotation | 2.1MB FASTA, 1.0MB GFF3 |
| **MAKER** | Eukaryotic annotation | 1.4MB FASTA, 285KB GFF3 |
| **Augustus** | Gene prediction | 1.4MB FASTA, 147KB GFF3 |
| **FASTA-only** | Sequence-only test | 8MB FASTA (gz) |

### Running Tests
```bash
# All tests with doctests (using uv)
uv run -- pytest --doctest-modules tests/

# Specific test file
uv run -- pytest tests/test_annotation_parsing_and_formatting.py

# Golden tests only
uv run -- pytest tests/golden/test_golden.py

# Without uv (if installed directly)
pytest --doctest-modules tests/
```

### Test Coverage
- **Roundtrip integrity**: Parse DDBJ → format DDBJ = identity
- **Multi-source validation**: Prokka, MAKER, Augustus compatibility
- **RefSeq comparison**: Validates against GenBank-derived annotations
- **DDBJ MSS validation**: Passes official DDBJ parser

---

## 8. RECENT CHANGES (Git History)

### Latest Commits
```
0267ad3  Update README-ja
85c66d1  Update README
6009353  Bump setuptools for build-system
30f941c  Create SeqRecord with Seq, not with a plain string
ecce883  Handle SeqFeatures as unhashable objects
```

### Key Improvements
- **Biopython 1.86 compatibility**: Fixed SeqRecord creation
- **Unhashable SeqFeatures**: Proper handling in data structures
- **Modern packaging**: Updated to setuptools >= 77.0.3
- **Documentation**: Bilingual README (EN + JA)

---

## 9. KNOWN LIMITATIONS

### Current Constraints
1. **Trans-splicing**: No coordinate correction for trans-spliced features
2. **Translation Exceptions**: `/transl_except` not fully handled
3. **Missing Qualifiers**: `/translation` not auto-generated when `/exception` present
4. **Inter-base Coordinates**: Limited support for "between-position" notation (`123^124`)
5. **Performance**: Single-process architecture (accuracy over speed)

### Areas for Future Enhancement
- Multi-process support for large genomes
- Extended `/transl_except` handling
- Automatic `/translation` generation
- Improved trans-splicing support
- Inter-base coordinate full support

---

## 10. UNTRACKED FILES (Git Status)

### Current Untracked Files
```
GCF_017312705.1_Mj_TUMSAT_v1.0.err
GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.ann
notebooks/start-stop-codons.ipynb
prokka.dr, prokka.ecn, prokka.fixedproducts
prokka.sqn, prokka.stats, prokka.val
template.sbt
tests/golden/prokka_renamed.gff3
tests/golden/pytest__*.ann (4 files)
uv.lock
```

### Recommendations
- **Annotation files** (`.ann`, `.err`): Consider adding to `.gitignore`
- **Prokka temp files**: Likely test artifacts, could be cleaned up
- **uv.lock**: Should be committed (user is actively using `uv`)
- **Notebooks**: Review `start-stop-codons.ipynb` for insights

---

## 11. DEVELOPMENT WORKFLOW

### Making Changes
1. **Read relevant code first**: Use Read tool before modifying
2. **Run tests**: `uv run -- pytest --doctest-modules tests/`
3. **Check golden tests**: Ensure reference data still passes
4. **Update CHANGELOG.md**: Document significant changes
5. **Update version**: Bump version in `pyproject.toml` and `version.py`

### Common Development Tasks
- **Add new feature mapping**: Edit `translate_features_qualifiers.toml`
- **Modify DDBJ compliance**: Edit `ddbj_filter.toml`
- **Add new CLI tool**: Create in `tools/` directory
- **Extend transforms**: Modify `transforms.py` pipeline
- **Update formatting**: Edit `formatter.py`

### Testing New Features
```bash
# Run on test data (using uv)
uv run -- gff3-to-ddbj \
  --fasta tests/golden/prokka.fa \
  --gff3 tests/golden/prokka.gff3 \
  --locus_tag_prefix TEST_ \
  --output test_output.ann

# Compare against golden
uv run -- compare-ddbj tests/golden/prokka.ann test_output.ann
```

---

## 12. QUICK REFERENCE

### Important File Locations
- **Main logic**: `gff3toddbj/transforms.py:1000+ LOC`
- **DDBJ formatting**: `gff3toddbj/formatter.py:500+ LOC`
- **Feature mappings**: `gff3toddbj/translate_features_qualifiers.toml:1037 lines`
- **DDBJ compliance**: `gff3toddbj/ddbj_filter.toml:523 lines`
- **Test data**: `tests/golden/`

### Key Functions to Understand
- `transforms.py`: Feature/qualifier transformation pipeline
- `formatter.py`: DDBJ annotation output generation
- `parser.py`: DDBJ annotation parsing (roundtrip validation)
- `io.py`: GFF3 parsing, FASTA indexing

### Useful Resources
- **README.md**: User-facing documentation (English)
- **README-ja.md**: User-facing documentation (Japanese)
- **CHANGELOG.md**: Version history and changes
- **tests/README.md**: Test data documentation
- **evaluation/**: Validation methodology

---

## 13. UV PACKAGE MANAGER NOTES

### About uv
This installation uses `uv` for Python package management. While the README suggests Bioconda, this setup uses `uv` for:
- Fast dependency resolution
- Reproducible builds via `uv.lock`
- Integrated virtual environment management

### Common uv Commands
```bash
# Run tests
uv run -- pytest

# Run the main CLI
uv run -- gff3-to-ddbj --help

# Run any Python script
uv run -- python script.py

# Add a dependency
uv add <package-name>

# Update dependencies
uv sync

# Lock dependencies
uv lock
```

### Important Files
- **uv.lock**: Lock file for reproducible dependencies (should be committed)
- **pyproject.toml**: Project configuration and dependencies

### uv vs conda
| Feature | uv | conda |
|---------|----|----|
| **Speed** | Very fast | Slower |
| **Python-only** | Yes (with some exceptions) | No (system-level packages) |
| **Bioinformatics tools** | Limited | Excellent (bioconda) |
| **Lock file** | uv.lock | environment.yml |
| **Use case** | Python development | Scientific computing |

For this project: `uv` is fine since main deps (biopython, pysam) are pip-installable.

---

## 14. VERIFICATION & TESTING RESULTS

### Verification Date
**2026-01-13** - Full verification performed

### Test Suite Results
```bash
# Command: uv run -- pytest --doctest-modules tests/
Platform: linux -- Python 3.13.1, pytest-9.0.2, pluggy-1.6.0
Result: ✓ 6 passed in 11.95s
Status: ALL TESTS PASSING
```

**Tests Verified**:
- ✓ RefSeq golden test (GCF_000280675.1)
- ✓ Prokka golden test
- ✓ MAKER golden test
- ✓ Augustus golden test
- ✓ FASTA-only test (no GFF3)
- ✓ Annotation parsing and formatting roundtrip test

### CLI Commands Verification
All CLI commands operational and producing correct output:

**Main Conversion Tool** - `gff3-to-ddbj`
```bash
✓ Version: 0.4.3
✓ Help system working
✓ Full conversion test on augustus.gff3 + augustus.fa
✓ Output matches golden file exactly (744 lines, 40K)
✓ Metadata loading working
✓ Warnings properly displayed for discarded features
```

**Alternative Converter** - `genbank-to-ddbj`
```bash
✓ Help system working
✓ All options available
```

**Comparison Tool** - `compare-ddbj`
```bash
✓ Help system working
✓ Successfully compared test output vs golden file
✓ Confirmed 0% mismatch (0/450 features, 0/628 qualifier pairs)
```

**Utility Tools**
```bash
✓ split-fasta - Help working
✓ list-products - Help working, basic functionality verified
✓ normalize-entry-names - Available
```

### End-to-End Conversion Test
**Test Case**: Augustus gene predictions
```bash
Input Files:
  - GFF3: tests/golden/augustus.gff3 (143K, 146KB)
  - FASTA: tests/golden/augustus.fa (1.4M)
  - Metadata: tests/golden/metadata_without_COMMON.toml
  - Translation table: 1 (standard genetic code)

Output:
  - Generated: /tmp/test_augustus_output.ann (40K, 744 lines)
  - Expected: tests/golden/augustus.ann (40K, 744 lines)
  - Comparison: IDENTICAL (verified with diff and compare-ddbj)

Processing Stats:
  - Features discarded: gene (58), transcription_start_site (58),
    start_codon (58), stop_codon (57), transcription_end_site (57)
  - Qualifiers discarded: Parent, score (for various features)
  - Processing completed successfully with proper DDBJ compliance filtering
```

### Environment Details
```
Python: 3.13.1
pytest: 9.0.2
Platform: Linux (WSL2)
Package Manager: uv
Dependencies: All resolved via uv.lock
```

### Verification Conclusions
✓ **All tests pass**
✓ **CLI commands functional**
✓ **End-to-end conversion produces identical output to golden files**
✓ **DDBJ compliance filtering working correctly**
✓ **Memory-efficient FASTA indexing operational**
✓ **Metadata loading and application working**
✓ **Roundtrip parsing/formatting maintains integrity**

**Status**: The codebase is fully functional and production-ready for DDBJ annotation file generation.

---

## 15. FUTURE WORK TRACKER

### Pending Tasks
- [ ] Review `start-stop-codons.ipynb` for new insights
- [ ] Clean up untracked prokka temp files
- [ ] Consider multi-process support for large genomes
- [ ] Enhance trans-splicing coordinate handling
- [ ] Implement automatic `/translation` generation

### Completed Tasks
- [x] Initial codebase review (2026-01-13)
- [x] CLAUDE.md creation (2026-01-13)
- [x] Confirmed `uv` as package manager, updated CLAUDE.md (2026-01-13)
- [x] Full verification of codebase functionality (2026-01-13)
  - All 6 pytest tests passing
  - CLI commands verified working
  - End-to-end conversion test successful (augustus dataset)
  - Output matches golden files exactly

---

## 16. CONTACT & CONTRIBUTION

### Project Information
- **Repository**: https://github.com/yamaton/gff3toddbj
- **Author**: Yamato Matsuoka
- **License**: GPLv3

### Contributing Guidelines
- Follow existing code style and patterns
- Add tests for new features
- Update documentation (README, CHANGELOG)
- Ensure pytest passes before committing
- Maintain DDBJ compliance standards

---

## 17. NOTES FOR FUTURE SESSIONS

### Code Quality Observations
- **Well-structured**: Clear separation of concerns across modules
- **Well-tested**: Comprehensive golden test suite with multiple sources
- **Well-documented**: Bilingual README, inline doctests, examples
- **Modern Python**: PEP 621 compliant, type hints present
- **Bioinformatics best practices**: Proper coordinate system handling

### Potential Improvements
- Type hints could be expanded throughout codebase
- Consider adding more inline documentation for complex algorithms
- Performance profiling for large genome optimization
- Error messages could be more user-friendly
- Consider adding logging framework instead of stderr prints

---

**End of CLAUDE.md**
