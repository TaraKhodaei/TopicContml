#!/usr/bin/env python
# coding: utf-8
import os
import sys
import shutil

def myparser():
    import argparse
    parser = argparse.ArgumentParser(description='convert fasta files to TopicContml loci file formaats')
    parser.add_argument('-f','--fasta_folder', dest='fasta_folder',
                        action='store', default='fasta', type=str,
                        help='the folder that contains the data (fasta in separate files). ')

                        
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args
    
    
    
    
#----------------------------------
def parse_fasta_file(input_file):
    """Return a dict of {id:gene_seq} pairs based on the sequences in the input FASTA file
    input_file -- a file handle for an input fasta file
    """
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []

    for line in input_file:
        line = line.strip()

        if line.startswith(">"):
            if curr_seq_id is not None:
                parsed_seqs[curr_seq_id] = ''.join(curr_seq)

            curr_seq_id = line[1:]
            curr_seq = []
            continue

        curr_seq.append(line)

    #Add the final sequence to the dict
    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
    return parsed_seqs



#----------------------------------
def output_fasta_file(parsed_seqs, fasta_clean_path, locus_num):
    """
    Transform FASTA file and write the output to a file with each line having an id,
    and then the sequence with 20 tab characters from begining
    """
    output_filepath = os.path.join(fasta_clean_path, locus_num)
    output_file = open(output_filepath,'w+')
    
    seqs_num = len(parsed_seqs)
    seqs_len = len(list(parsed_seqs.values())[0])
    

    output_file.write('{} {}\n'.format(seqs_num, seqs_len))

    for seq_id,seq in parsed_seqs.items():
        seq_id_name = seq_id.split()[1:3]
        seq_name = str(seq_id_name[0][0])+'.'+ str(seq_id_name[1])
        output_file.write(f'{seq_name:<20}{"".join(map(str, seq))}\n')
#    return output_file



#----------------------------------
if __name__ == '__main__':
    args = myparser() # parses the commandline arguments
    fasta_folder = args.fasta_folder
    fasta_folder_clean = 'loci_'+fasta_folder

    parent_dir = os.getcwd()
    fasta_path = os.path.join(parent_dir,fasta_folder)    #path of fasta_folder


    #creating [fasta folder name]_clean
    fasta_clean_path = os.path.join(parent_dir, fasta_folder_clean)
    dir = fasta_folder_clean
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.mkdir(fasta_clean_path)

    #os.listdir() method in python is used to get the list of all files and directories in the specified directory.
    for num,file in enumerate(os.listdir(fasta_path)):
        locus_num = 'locus'+str(num)+'.txt'
        fasta_file = os.path.join(fasta_path,file)
        print(f"{file} ---> {locus_num}")
        with open(fasta_file, 'r' ) as myfile:
            parsed_seqs = parse_fasta_file(myfile)
        output_fasta_file(parsed_seqs, fasta_clean_path, locus_num)
    
    
 

