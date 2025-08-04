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
def find_orfs(sequence, min_len, strand="+"):
    """Strand-specific ORF detection that allows adjacent, non-overlapping ORFs and blocks true overlaps."""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    results = []
    seq_len = len(sequence)

    used = [False] * seq_len  # Used bases for this strand only

    for frame in range(3):
        pos = frame
        while pos <= seq_len - 3:
            codon = sequence[pos:pos+3]

            if codon == start_codon:
                # Check if this start codon overlaps an existing ORF
                if any(used[pos:pos+3]):
                    pos += 1
                    continue

                end = pos + 3
                while end <= seq_len - 3:
                    stop = sequence[end:end+3]
                    if stop in stop_codons:
                        orf_len = end + 3 - pos
                        if orf_len >= min_len * 3:
                            # Only block ORFs that truly overlap
                            if not any(used[pos:end]):
                                orf_seq = sequence[pos:end+3]
                                if strand == '+':
                                    start_pos = pos + 1
                                    frame_num = frame + 1
                                else:
                                    start_pos = seq_len - pos
                                    frame_num = frame + 4
                                results.append((frame_num, start_pos, orf_seq, strand))
                                for i in range(pos, end):  # exclude stop codon
                                    used[i] = True
                                pos = end + 3  # move past stop codon
                                break
                    end += 3
                else:
                    pos += 3
            else:
                pos += 1
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
    # Ask user for FASTA file and minimum ORF length
    fasta_file = input("Enter FASTA filename: ")
    min_len_input = input("Enter minimum ORF length (default 5): ")
    min_len = int(min_len_input) if min_len_input.isdigit() else 5

    # Load sequences from the FASTA file
    sequences = load_fasta(fasta_file)
    print(f"Loaded sequences: {list(sequences.keys())}")

    # Set paths for output FASTA file and visualization images
    output_path = "./output/orfs/orf_output.fasta"
    visualization_path = "./output/visualization/orf_visualization.png"

    # Make sure output directories exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.makedirs(os.path.dirname(visualization_path), exist_ok=True)

    # Store all ORFs across all sequences for summary and visualization
    all_orfs = []

    # Open output file for writing ORF results
    with open(output_path, 'w') as out_file:
        # Loop through each sequence in the FASTA file
        for header, seq in sequences.items():
            # Find forward strand ORFs
            orfs_forward = find_orfs(seq, min_len, strand='+')
            all_orfs.extend([(header, *orf) for orf in orfs_forward])

            # Find reverse strand ORFs (from reverse complement)
            rev_seq = reverse_complement(seq)
            orfs_reverse = find_orfs(rev_seq, min_len, strand='-')
            all_orfs.extend([(header, *orf) for orf in orfs_reverse])

            # Print and write formatted ORF output
            for orf in orfs_forward + orfs_reverse:
                frame, position, orf_seq, strand = orf
                formatted = format_orf_output(header, frame, position, orf_seq, strand)
                print(formatted.strip())
                out_file.write(formatted)

    # Create and save visualization charts
    create_visualization(all_orfs, visualization_path)

# Run the program
if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
