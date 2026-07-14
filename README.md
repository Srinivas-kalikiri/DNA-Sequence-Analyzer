# DNA Sequence Analyzer

A command-line Python tool for analyzing DNA sequences — nucleotide composition, GC content, complement/reverse complement, transcription, translation, and motif search. Supports both raw sequence input and FASTA files.

## Features
- **Validation** — rejects sequences containing anything other than A, T, G, C with a clear error message
- **Composition analysis** — nucleotide frequency counts and GC content percentage
- **Complement & reverse complement** generation
- **Transcription** — DNA → RNA (T → U)
- **Translation** — DNA → protein using the standard genetic code, stopping at the first stop codon
- **Reading frames** — translates all 3 forward and 3 reverse reading frames
- **Motif search** — finds every position of a subsequence/motif within a larger sequence
- **FASTA file support** — parses multi-record FASTA files and analyzes each sequence in turn

## Tech Stack
- Python 3.9+ (standard library only — no external dependencies)

## Project Structure
```
dna-sequence-analyzer/
├── dna_analyzer.py           # Core library + CLI
├── sample_data/
│   └── sample.fasta          # Example FASTA file for testing
├── tests/
│   └── test_dna_analyzer.py  # Unit tests
├── requirements.txt
└── README.md
```

## Setup
```bash
git clone https://github.com/<Srinivas-kalikiri>/dna-sequence-analyzer.git
cd dna-sequence-analyzer
```
No external dependencies are required — just Python 3.9+.

## Usage

Analyze a raw sequence:
```bash
python dna_analyzer.py --seq ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

Analyze every sequence in a FASTA file:
```bash
python dna_analyzer.py --fasta sample_data/sample.fasta
```

Search for a motif:
```bash
python dna_analyzer.py --seq ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG --motif GGGCC
```

### Example Output
```
=== Sequence ===
Length: 39 bp
Nucleotide frequency: {'A': 9, 'T': 8, 'G': 14, 'C': 8}
GC content: 56.41%
Complement:         TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
RNA transcript:     AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein (frame 1):  MAIVMGR
Motif 'GGGCC' found at position(s): [14]
```

## How It Works
The genetic code lookup uses the standard 64-codon table mapping each codon to its one-letter amino acid code (or a stop symbol `*`). Translation reads the sequence three bases at a time starting from the first base and stops as soon as a stop codon is encountered, mirroring how ribosomes terminate translation in real biology.

## Running Tests
```bash
python -m unittest discover tests
```

## Possible Improvements
- Support ambiguous IUPAC nucleotide codes (N, R, Y, etc.)
- ORF (Open Reading Frame) finder across all 6 reading frames
- Restriction enzyme cut-site detection
- Export analysis results to CSV/JSON

