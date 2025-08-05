from Bio.Seq import Seq
import matplotlib.pyplot as plt
import os

def load_fasta(filepath):
    # Team Member Name: Harrison Foncha
    # Parse a multi-line FASTA file, return dictionary {header: sequence}
    sequences = {}
    with open(filepath, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = ''.join(seq_lines)
                header = line[1:]  # remove ">"
                seq_lines = []
            else:
                seq_lines.append(line.upper())
        if header:
            sequences[header] = ''.join(seq_lines)
    return sequences

def reverse_complement(seq):
    # Team Member Name: Harrison
    # Return reverse complement of sequence using Bio.Seq
    return str(Seq(seq).reverse_complement())

def find_orfs(header, sequence, min_len=100, strand="+"):
    # Team Member Name: Harrison
    # Identify ORFs in all 3 reading frames
    stop_codons = ['TAA', 'TAG', 'TGA']
    start_codon = 'ATG'
    orfs = []

    for frame in range(3):
        start = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == start_codon and start is None:
                start = i
            elif codon in stop_codons and start is not None:
                orf_seq = sequence[start:i+3]
                if len(orf_seq) >= min_len:
                    orfs.append({
                        "header": header,
                        "frame": f"{strand}{frame+1}",
                        "start": start + 1,
                        "end": i + 3,
                        "length": len(orf_seq),
                        "sequence": orf_seq
                    })
                start = None
    return orfs

def format_orf_output(header, frame, position, seq):
    # Team Member Name: Harrison
    # Return formatted FASTA header and codon-separated sequence
    codons = ' '.join([seq[i:i+3] for i in range(0, len(seq), 3)])
    new_header = f">{header} | Frame: {frame} | Start: {position}"
    return f"{new_header}\n{codons}"

def create_visualization(orf_data, output_path):
    # Team Member Name: Harrison
    # Create a visualization of ORF lengths and save to file
    lengths = [orf["length"] for orf in orf_data]
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=20, color='skyblue', edgecolor='black')
    plt.title("Distribution of ORF Lengths")
    plt.xlabel("ORF Length (bp)")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    # Team Member Name: Harrison
    filepath = r"C:\Users\Harrison West\Assignment\sequence.fasta"  # Update filename if needed
    output_orf_file = r"C:\Users\Harrison West\Assignment\orf_output.fasta"
    output_plot = r"C:\Users\Harrison West\Assignment\orf_plot.png"
    min_orf_len = 100

    sequences = load_fasta(filepath)
    all_orfs = []

    with open(output_orf_file, 'w') as out:
        for header, seq in sequences.items():
            # Forward strand
            orfs_fwd = find_orfs(header, seq, min_orf_len, "+")
            # Reverse strand
            rc_seq = reverse_complement(seq)
            orfs_rev = find_orfs(header, rc_seq, min_orf_len, "-")
            all_orfs.extend(orfs_fwd + orfs_rev)

            for orf in orfs_fwd + orfs_rev:
                out.write(format_orf_output(
                    orf["header"], orf["frame"], orf["start"], orf["sequence"]
                ) + "\n")

    create_visualization(all_orfs, output_plot)
    print(f"✔ ORF analysis complete. Results saved to:\n- {output_orf_file}\n- {output_plot}")

if __name__ == "__main__":
    main()