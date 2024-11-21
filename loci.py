import os
import sys
import math

# Function to format sequences in single line
def format_fasta_single_line(input_file):
    with open(input_file, 'r') as file:
        data = file.read()

    fasta_entries = data.strip().split(">")
    fasta_formatted = []
    
    for entry in fasta_entries:
        if entry.strip():  # Ensure it's not empty
            entry_parts = entry.split("\n", 1)  # Split by first newline
            name = entry_parts[0].strip()
            seq = entry_parts[1].replace("\n", "").strip()  # Remove all newlines from the sequence
            fasta_formatted.append(f">{name}\n{seq}")  # Keep sequence in one line
    
    formatted_fasta = "\n".join(fasta_formatted)

    # Overwrite the input file with the formatted content
    with open(input_file, 'w') as file:
        file.write(formatted_fasta)

def myparser():
    import argparse
    parser = argparse.ArgumentParser(description='Generate loci from Fasta files')
    
    parser.add_argument('-fr','--fasta_unassebbled_reads', dest='fasta_unassebbled_reads',
                        default=None, action='store',type=int,
                        help='Specifies the number of loci (n) to generate from unassembled FASTA read files. The value of n is provided directly with the -fr option (e.g., -fr 5 for 5 loci). If this option is used, the script treats the FASTA files as unassembled reads. For each file, it concatenates every ceil(total_reads_per_file / n) reads, then generates n loci files by combining the concatenated sequences from all FASTA files.')
    parser.add_argument('-fg','--fasta_complete_genome', dest='fasta_complete_genome',
                        default=None, action='store',type=int,
                        help='Specifies that the input FASTA files are complete genomes, with each file containing a single genome sequence. If this option is used, the script generates n loci by dividing the sequence in each FASTA file into n segments. It then creates n loci files by combining corresponding segments from all FASTA files.')
    parser.add_argument('-if','--input_folder', dest='input_folder_path',
                        default=None, action='store',  type=str,
                        help='Specifies the name of the input folder containing the FASTA files to be used.')
    parser.add_argument('-of','--output_folder', dest='output_folder_path',
                        default=None, action='store',type=str,
                        help='Specifies the name of the output folder where the generated loci files will be saved.')
                        
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    return args

def read_sequences(file_path):
    with open(file_path, 'r') as file:
        while True:
            lines = [file.readline().strip() for _ in range(2)]
            if not lines[0]:
                break
            yield lines[1]

def concatenate_sequences(sequences, chunk_size):
    chunk = []
    for seq in sequences:
        chunk.append(seq)
        if len(chunk) == chunk_size:
            yield ''.join(chunk)
            chunk = []
    if chunk:
        yield ''.join(chunk)

def format_sequences(sequences, seq_names):
    formatted_sequences = []
    for idx, seq in enumerate(sequences):
        prefix = f"{seq_names[idx]:<11}"      #phylip standard format
        formatted_seq = f"{prefix[:10]} {seq}"
        formatted_sequences.append(formatted_seq)
    return formatted_sequences

def write_sequences(filename, header, sequences):
    with open(filename, 'w') as file:
        file.write(header + '\n')
        for seq in sequences:
            file.write(seq + '\n')


#====================================================================================
if __name__ == "__main__":

    args = myparser() # parses the commandline arguments
    fasta_unassebbled_reads = args.fasta_unassebbled_reads
    fasta_complete_genome = args.fasta_complete_genome
    input_folder_path = args.input_folder_path
    output_folder_path = args.output_folder_path

    input_files_names = []
    
    # Sort and filter input files
    for file_name in sorted(os.listdir(input_folder_path)):
        # Ignore hidden files like .DS_Store and non-FASTA files
        if not file_name.startswith('.') and os.path.isfile(os.path.join(input_folder_path, file_name)):
            input_files_names.append(file_name)
    
    print(f"\ninput_files_names = {input_files_names}")

    files_dict = {}
    
    # Step to format each FASTA file to ensure sequences are on a single line
    for file_name in input_files_names:
        file_path = os.path.join(input_folder_path, file_name)
        format_fasta_single_line(file_path)  # Ensure sequences in the file are single-line
    
    if fasta_unassebbled_reads:
        num_loci = fasta_unassebbled_reads
        for idx, file in enumerate(input_files_names):
            file_path = os.path.join(input_folder_path, file)
            
            # Convert generator to list once to get all sequences
            sequences = list(read_sequences(file_path))

            # Calculate chunk size, ensuring it is at least 1
            chunk_size = max(1, math.floor(len(sequences) / num_loci))

            # Concatenate sequences using the calculated chunk size
            concatenated_sequences = list(concatenate_sequences(sequences, chunk_size))
            # Save the sequences in a dictionary
            files_dict[idx] = concatenated_sequences
        
    if fasta_complete_genome:
        num_loci = fasta_complete_genome
        for idx, file in enumerate(input_files_names):
            file_path = os.path.join(input_folder_path, file)
            
            # Convert generator to list once to get all sequences
            sequences = list(read_sequences(file_path))

            # Check if the file contains at least one sequence
            if len(sequences) > 0:
                full_sequence = sequences[0]  # Get the complete genome sequence from the file

                # Calculate chopping size, ensuring it is at least 1
                chop_size = max(1, math.floor(len(full_sequence) / num_loci))

                # Generate exactly `num_loci` segments
                chopped_gene = [full_sequence[i * chop_size: (i + 1) * chop_size] for i in range(num_loci - 1)]
                chopped_gene.append(full_sequence[(num_loci - 1) * chop_size:])

                # Save the chopped gene in a dictionary
                files_dict[idx] = chopped_gene
            else:
                print(f"No sequences found in {file_path}.")
        
                
        
    # Create output files by combining files_dict sequences from each input file
    output_dir = output_folder_path
    os.makedirs(output_dir, exist_ok=True)
    
    for i in range(num_loci):
        combined_sequences = [files_dict[file_idx][i] for file_idx in range(len(input_files_names))]
        formatted_sequences = format_sequences(combined_sequences, input_files_names)
        
        num_sequences = len(formatted_sequences)
        lengths = "|".join(str(len(seq)) for seq in combined_sequences)
        header = f"{num_sequences} {lengths}"
        
        output_filename = os.path.join(output_dir, f'locus{i}.txt')
        write_sequences(output_filename, header, formatted_sequences)
        
    print(f"\nProcessing complete. Output loci created in folder '{output_dir}'.")
