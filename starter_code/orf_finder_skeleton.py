# BIFS617 ORF Finder Skeleton
# Team members:
# Date:
# Fill in your code where prompted.

def load_fasta(filepath):
    # Team Member Name:
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    pass

def reverse_complement(seq):
    # Team Member Name:
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)
    pass

def find_orfs(header, sequence, min_len, strand="+"):
    # Team Member Name:    
    # TODO: Identify ORFs in all 3 reading frames for one strand
    pass

def format_orf_output(header, frame, position, seq):
    # Team Member Name:
    # TODO: Return formatted FASTA header and codon-separated sequence
    pass

def create_visualization(orf_data, output_path):
    # Team Member Name:
    # TODO: create a visualization, save the file, for your ORF output
    pass

def main():
    # Team Member Name:
    # TODO: Implement user input, sequence processing, and ORF printing and save the file
    pass

if __name__ == "__main__":
    main()
