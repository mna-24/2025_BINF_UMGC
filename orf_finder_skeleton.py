# BIFS617 ORF Finder Skeleton
# Team members: Raquel, 
# Date:
# Fill in your code where prompted.

import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict

# Team Member Name: Raquel

# TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
def load_fasta(filepath):
    """Load sequences from a multi-FASTA file into a dictionary."""
    sequences = {}
    header = None
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:].strip()
                sequences[header] = ''
            elif header:
                sequences[header] += line.upper()
    return sequences

# TODO: Return reverse complement of sequence (optional: use Bio.Seq)
def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ATCGN', 'TAGCN')
    return seq.translate(complement)[::-1]

# TODO: Identify ORFs in all 3 reading frames for one strand
def find_orfs(header, sequence, min_len, strand="+"):
    """Identify ORFs in all 3 reading frames of one strand."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    results = []

    for frame in range(3):
        pos = frame
        while pos <= len(sequence) - 3:
            codon = sequence[pos:pos+3]
            if codon == start_codon:
                for end in range(pos + 3, len(sequence), 3):
                    stop = sequence[end:end+3]
                    if stop in stop_codons:
                        orf_seq = sequence[pos:end+3]
                        if len(orf_seq) >= min_len * 3:
                            if strand == '-':
                                start_pos = len(sequence) - pos
                            else:
                                start_pos = pos + 1
                            frame_num = frame + 1 if strand == '+' else frame + 4
                            results.append((frame_num, start_pos, orf_seq, strand))
                        break  # Do NOT skip ahead — allow overlapping ORFs
                pos += 1  # Slide one base to allow overlapping ORFs
            else:
                pos += 3
    return results

# TODO: Return formatted FASTA header and codon-separated sequence
def format_orf_output(header, frame, position, seq, strand):
    direction = 'FOR' if strand == '+' else 'REV'
    formatted = f">{header} | FRAME = {frame} | POS = {position} | LEN = {len(seq)} | {direction}\n"
    codon_chunks = ' '.join([seq[i:i+3] for i in range(0, len(seq), 3)])
    return formatted + codon_chunks + "\n"

# TODO: create a visualization, save the file, for your ORF output
def create_visualization(orf_data, output_path):
    """Create a bar chart showing number of ORFs per sequence per strand and a histogram of ORF lengths."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    summary = defaultdict(lambda: {'FOR': 0, 'REV': 0})
    orf_lengths = []

    for entry in orf_data:
        header, frame, position, seq, strand = entry
        direction = 'FOR' if strand == '+' else 'REV'
        summary[header][direction] += 1
        orf_lengths.append(len(seq))

    labels = list(summary.keys())
    forward_counts = [summary[label]['FOR'] for label in labels]
    reverse_counts = [summary[label]['REV'] for label in labels]

    x = range(len(labels))
    plt.figure(figsize=(10, 5))
    plt.bar(x, forward_counts, width=0.4, label='Forward', align='center')
    plt.bar(x, reverse_counts, width=0.4, label='Reverse', align='edge')
    plt.xlabel('Sequence ID')
    plt.ylabel('ORF Count')
    plt.title('Number of ORFs per Sequence by Strand')
    plt.xticks(x, labels, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(output_path), "orf_count_by_strand.png"))
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.hist(orf_lengths, bins=10, color='skyblue', edgecolor='black')
    plt.xlabel('ORF Length')
    plt.ylabel('Frequency')
    plt.title('Distribution of ORF Lengths')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# TODO: Implement user input, sequence processing, and ORF printing and save the file
def main():
    fasta_file = input("Enter FASTA filename: ")
    min_len_input = input("Enter minimum ORF length (default 5): ")
    min_len = int(min_len_input) if min_len_input.isdigit() else 5

    sequences = load_fasta(fasta_file)
    print(f"Loaded sequences: {list(sequences.keys())}")

    output_path = "./output/orfs/orf_output.fasta"
    visualization_path = "./output/visualization/orf_visualization.png"

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.makedirs(os.path.dirname(visualization_path), exist_ok=True)

    all_orfs = []
    with open(output_path, 'w') as out_file:
        for header, seq in sequences.items():
            orfs_forward = find_orfs(header, seq, min_len, strand='+')
            all_orfs.extend([(header, *orf) for orf in orfs_forward])

            rev_seq = reverse_complement(seq)
            orfs_reverse = find_orfs(header, rev_seq, min_len, strand='-')
            all_orfs.extend([(header, *orf) for orf in orfs_reverse])

            for orf in orfs_forward + orfs_reverse:
                frame, position, orf_seq, strand = orf
                formatted = format_orf_output(header, frame, position, orf_seq, strand)
                print(formatted.strip())
                out_file.write(formatted)

    create_visualization(all_orfs, visualization_path)

if __name__ == "__main__":
    main()
