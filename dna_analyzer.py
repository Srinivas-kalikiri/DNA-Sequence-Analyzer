"""
DNA Sequence Analyzer
A command-line tool for analyzing DNA sequences: composition, GC content,
complement/reverse complement, transcription, translation, and motif search.
"""

import argparse
import re
from collections import Counter

VALID_BASES = set("ATGC")

COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Standard genetic code: codon -> single-letter amino acid ('*' = stop)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


class InvalidSequenceError(ValueError):
    pass


def clean_sequence(seq):
    return re.sub(r"\s", "", seq).upper()


def validate_sequence(seq):
    if not seq:
        raise InvalidSequenceError("Sequence is empty.")
    invalid_chars = sorted(set(seq) - VALID_BASES)
    if invalid_chars:
        raise InvalidSequenceError(
            f"Sequence contains invalid character(s): {', '.join(invalid_chars)}. "
            f"Only A, T, G, C are allowed."
        )
    return True


def nucleotide_frequency(seq):
    seq = clean_sequence(seq)
    validate_sequence(seq)
    counts = Counter(seq)
    return {base: counts.get(base, 0) for base in "ATGC"}


def gc_content(seq):
    seq = clean_sequence(seq)
    validate_sequence(seq)
    gc = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


def complement(seq):
    seq = clean_sequence(seq)
    validate_sequence(seq)
    return "".join(COMPLEMENT_MAP[base] for base in seq)


def reverse_complement(seq):
    return complement(seq)[::-1]


def transcribe(seq):
    seq = clean_sequence(seq)
    validate_sequence(seq)
    return seq.replace("T", "U")


def translate(seq):
    """Translates a DNA sequence into protein, stopping at the first stop codon."""
    seq = clean_sequence(seq)
    validate_sequence(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "*":
            break
        protein.append(amino_acid)
    return "".join(protein)


def find_motif(seq, motif):
    """Returns every 0-indexed start position where motif occurs in seq (overlaps included)."""
    seq = clean_sequence(seq)
    motif = clean_sequence(motif)
    validate_sequence(seq)
    validate_sequence(motif)
    positions = []
    start = 0
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def reading_frames(seq):
    """Translates all 3 forward and 3 reverse reading frames."""
    seq = clean_sequence(seq)
    validate_sequence(seq)
    frames = {}
    for i in range(3):
        frames[f"Frame +{i + 1}"] = translate(seq[i:])
    rc = reverse_complement(seq)
    for i in range(3):
        frames[f"Frame -{i + 1}"] = translate(rc[i:])
    return frames


def parse_fasta(filepath):
    """Parses a FASTA file into a dict of {header: sequence}."""
    sequences = {}
    header = None
    seq_chunks = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            sequences[header] = "".join(seq_chunks)
    return sequences


def analyze_sequence(seq, label="Sequence"):
    seq = clean_sequence(seq)
    print(f"\n=== {label} ===")
    print(f"Length: {len(seq)} bp")
    print(f"Nucleotide frequency: {nucleotide_frequency(seq)}")
    print(f"GC content: {gc_content(seq)}%")
    print(f"Complement:         {complement(seq)}")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print(f"RNA transcript:     {transcribe(seq)}")
    print(f"Protein (frame 1):  {translate(seq)}")


def main():
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a raw string or FASTA file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--seq", help="Raw DNA sequence, e.g. ATGCGTAC")
    group.add_argument("--fasta", help="Path to a FASTA file")
    parser.add_argument("--motif", help="Search for a motif/subsequence and report its positions")

    args = parser.parse_args()

    if args.seq:
        analyze_sequence(args.seq)
        if args.motif:
            positions = find_motif(args.seq, args.motif)
            print(f"Motif '{args.motif}' found at position(s): {positions if positions else 'not found'}")

    elif args.fasta:
        records = parse_fasta(args.fasta)
        for header, seq in records.items():
            analyze_sequence(seq, label=header)
            if args.motif:
                positions = find_motif(seq, args.motif)
                print(f"Motif '{args.motif}' found at position(s): {positions if positions else 'not found'}")


if __name__ == "__main__":
    main()
