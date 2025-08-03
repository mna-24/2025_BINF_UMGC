# BIFS617 ORF Finder Skeleton
# Team members:
# Date:
# Fill in your code where prompted.

def load_fasta(filepath):
    # Team Member Name:
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    with open(filepath, 'r') as file: #pass file path to function instead of putting it here, later!!
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

 #debug lines to see output
func = load_fasta('/Users/michaelalemu/PycharmProject/Group_project_BIFS_617/2025_BINF_UMGC/project_input/orfs_test_input.fasta') #calling the method and passing file path.
print(f"load fasta & parse > to dict function result: {func}")

from Bio.Seq import Seq
def reverse_complement(original_seq):
    # Team Member Name:
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)

    seq_dict = load_fasta(original_seq)
    for header,sequence in seq_dict.items():
        seq_obj = Seq(sequence)
        rev_seq = seq_obj.reverse_complement()
        rev_seq_strng = str(rev_seq)

        seq_dict[header] = rev_seq_strng

    return seq_dict

file_path = '/Users/michaelalemu/PycharmProject/Group_project_BIFS_617/2025_BINF_UMGC/project_input/orfs_test_input.fasta' #calling the method and passing file path.
func_2 = reverse_complement(file_path)
print(f"reverse complement function output: {func_2}")


def find_orfs():
    # Team Member Name:    
    # TODO: Identify ORFs in all 3 reading frames for one strand (The project objective says all 6 reading frames though???
    seq_dict = load_fasta('/Users/michaelalemu/PycharmProject/Group_project_BIFS_617/2025_BINF_UMGC/project_input/orfs_test_input.fasta') #calling the method and passing file path.

    # header, sequence, min_len, strand = "+"
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    # orfs = []
    orf_dict = {}

    # iterating through the sequence dictionary
    for header,sequence in seq_dict.items():

            orfs = [] #empty list to collect orfs. (ORF is nucleotide sequence that codes for protein.)

            for n in range(3): # for getting sequence of numbers 0 to 3, representing the reading frames, to be fed for the next loop. (here is where the shift happens)
                i = n
                while i < len(sequence) -2: # loops until we hit len of seq (-2) so that we don't fall out of the end of the seq.
                # find the start codon, grab its position/index
                    codon = sequence[i:i+3] #here we are reading in 3's, each codon (3 nucleotides) at a time,
                    if codon == start_codon: # find start codons (ATG) in this reading frame.
                        start_codon_index = i # save starting index of the start codon
                        i +=3 #move by 3, to look at the next codon after the start codon, in this reading frame

                        # finding the stop codon, grab its position/index
                        while i < len(sequence) -2: # i starts right after the start codon here (i += 3 from previous if statement)
                            codon = sequence[i:i+3] #Again here we are reading in 3's, each codon
                            if codon in stop_codons: # if we find a stop codon, grab the codon index.
                                stop_codon_index = i + 3 # since stop codons are included in the orf, we want the index of i + 3

                                orfs.append(sequence[start_codon_index:stop_codon_index]) #saving each ORF after we loop and find a start and end codon.
                                break # break the inner while loop and go back to the start codon while loop.

                            i +=3 # so that the next start codon search starts at the end of the stop codon in the seq.

                    else:
                        i += 3 # if the current codon is not a start codon, we want to go back to the outer while loop to find the next codon (so i is now i + 3).
            # return orfs  #DEBUG: At the end, this function returns the orfs in a list

            orf_dict[header] = orfs #we might move this along with the looping through the main seq_dict to the main method for reusability.
    return orf_dict

#debug lines
my_orf_finder_func = find_orfs()
print(f"orfs: {my_orf_finder_func}")




def format_orf_output():
    # Team Member Name:
    # TODO: Return formatted FASTA header and codon-separated (spaces between each codon) sequence
    #pos arguments: header, frame, position, seq
    orf_dict = find_orfs()
    for Header, Sequence in orf_dict.items():
        Header = f">{Header}"
        for i in Sequence:
            codons = []
            n = 0
            while n < len(i):
                codons.append(i[n:n+3])
                n += 3

            formatted_seq = ' '.join(codons)

    return Header + "\n" + formatted_seq





FASTA_format_orf= format_orf_output()
print(FASTA_format_orf)

def create_visualization(orf_data, output_path):
    # Team Member Name:
    # TODO: create a visualization, save the file, for your ORF output
    pass

def main():
    # Team Member Name:
    # TODO: Implement user input, sequence processing, and ORF printing and save the file
    # func = load_fasta('/Users/michaelalemu/PycharmProject/Group_project_BIFS_617/2025_BINF_UMGC/project_input/orfs_test_input.fasta') #calling the method and passing file path.

    #parse the seq dictionary, and pass header, sequence and pass it to def find_orfs function when I call it; so that I can also reuse it for the reverse complimentary sequence as well.
    pass

if __name__ == "__main__":
    main()


