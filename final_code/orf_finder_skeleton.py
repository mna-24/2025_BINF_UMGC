# BIFS617 ORF Finder Skeleton
# Team members: Michael, Raquel, Harrison
# Date: 05Aug25
# Fill in your code where prompted.
import os
from collections import defaultdict

from matplotlib import pyplot as plt


def load_fasta(filepath):
    # Team Member Name:
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    with open(filepath, 'r') as file: #pass file path to function
        seq_dict = {} #initialize empty dictionary
        header = None # set my variable (header) to None as placeholder
        seq = [] #init empty list since the seq is multiline, save each iteration and will join it later before saving to the dict.
        for line in file: #looping through the content of the file, line-by-line by iterating.
            line = line.strip() #remove any white space on new line space
            if line.startswith('>'): # each header starts with '>', which has its seq following it in subsequent lines
                if header: # if its a header line.
                    seq_dict[header] = ''.join(seq) #save header and seq from previous iteration
                header = line[1:] # setting the new header start after '>'
                seq = [] #reset seq to empty to prevent appending next seq on the previous one.
            else:
                line = line.upper()
                seq.append(line) # append this iterations sequence in seq list.

        if header: # saving the last header:seq pair from the loop since my save code is inside the loop, and the loop has ended when there are no more lines.
            seq_dict[header] = ''.join(seq)
    return seq_dict


from Bio.Seq import Seq
def reverse_complement(original_seq):
    # Team Member Name:
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)

    seq_obj = Seq(original_seq) #creating Biopython seq object
    rev_seq = seq_obj.reverse_complement() #This generates the reverse compliment seq
    rev_seq_strng = str(rev_seq) #make the output a string


    return rev_seq_strng



def find_orfs(header, seq, min_len, strand = '+' ):
    # Team Member Name:
    # TODO: Identify ORFs in all 3 reading frames for one strand (The project objective says all 6 reading frames though???
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    results = [] #empty list to collect orfs. (ORF is nucleotide sequence that codes for protein.)

    for frame in range(3): #for the 3 reading frames
        pos = frame
        while pos <= len(seq) - 3: # loops until we hit len of seq (-2) so that we don't fall out of the end of the seq.
                # find the start codon, grab its position/index
            codon = seq[pos:pos+3] #here we are reading in 3's, each codon (3 nucleotides) at a time,
            if codon == start_codon: # find start codons (ATG) in this reading frame.
                for end in range(pos + 3, len(seq), 3): # looping until finding the stop codon
                    stop = seq[end:end+3]
                    if stop in stop_codons:
                        orf_seq = seq[pos:end+3]
                        if len(orf_seq) >= min_len:
                            #grab ORF start position
                            if strand == '-':
                                start_pos = len(seq) - pos
                            else:
                                start_pos = pos + 1
                            frame_num = frame + 1 if strand == '+' else frame + 4
                            results.append((header, frame_num, start_pos, orf_seq, strand))
                        break #stop after encountering first stop codon
                pos += 1 #This allows overlapping ORFs
            else:
                pos += 3 # move to the next codon
    return results




def format_orf_output(header, frame, position, seq, strand):
    # Team Member Name:
    # TODO: Return formatted FASTA header and codon-separated (spaces between each codon) sequence

    direction = 'FOR' if strand == '+' else 'REV'
    f_header = f"> {header} | Frame = {frame} | POS = {position} | LEN = {len(seq)} | {direction} " #concocate
    codons = ' '.join([seq[i:i+3] for i in range(0,len(seq), 3)])
    formatted_result = f"{f_header}\n{codons}\n"

    return formatted_result




def create_visualization(orf_data, output_path):
    # Team Member Name:
    # TODO: create a visualization, save the file, for your ORF output
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


 #Histogram of ORF lengths
    plt.figure(figsize=(8, 5))
    plt.hist(orf_lengths, bins=10, color='skyblue', edgecolor='black')
    plt.xlabel('ORF Length')
    plt.ylabel('Frequency')
    plt.title('Distribution of ORF Lengths')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    # Team Member Name:
    # TODO: Implement user input, sequence processing, and ORF printing and save the file

    fasta_file = input("Enter FASTA filename: ")
    min_len_input = input("Enter minimum ORF length (default 5): ")
    min_len = int(min_len_input) if min_len_input.isdigit() else 5 # if min_len not specified, default to 3

    #loading the sequence
    sequences = load_fasta(fasta_file)
    # print(f"Loaded sequences: {list(sequences.keys())}")

    #specifying output and visualization path
    output_path = "./output/orfs/orf_output.fasta"
    visualization_path = "./output/visualization/orf_visualization.png"

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.makedirs(os.path.dirname(visualization_path), exist_ok=True)


    all_orfs = []


    for header, seq in sequences.items():
        # Forward strand ORFs
        fwd_orfs = find_orfs(header, seq, min_len, strand='+')
        all_orfs.extend(fwd_orfs)

        # Reverse strand ORFs
        rev_seq = reverse_complement(seq)
        rev_orfs = find_orfs(header, rev_seq, min_len, strand='-')
        all_orfs.extend(rev_orfs)

    #write results (formatted ORFs) to the output file
    with open(output_path, 'w') as out_file:
        for orf in all_orfs:
            header, frame, position, orf_seq, strand = orf
            formatted = format_orf_output(header, frame, position, orf_seq, strand)
            print(formatted.strip())
            out_file.write(formatted)


    create_visualization(all_orfs, visualization_path)

    # ORF summary by sequence
    summary = defaultdict(lambda: {'FOR': 0, 'REV': 0})
    for header, frame, position, seq, strand in all_orfs:
        direction = 'FOR' if strand == '+' else 'REV'
        summary[header][direction] += 1

    # print("\nORF Summary:")
    for header in sequences:
        fwd = summary[header]['FOR']
        rev = summary[header]['REV']
        print(f"> {header}")
        print(f"{fwd} ORFs in forward, {rev} in reverse\n")


if __name__ == "__main__":
    main()


