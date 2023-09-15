#!/usr/bin/env python
# coding: utf-8

import phylip
import os
import sys
import shutil
import numpy as np

DEBUG=False


    
def myparser():
    import argparse
    parser = argparse.ArgumentParser(description='generates nexus file of concatenated loci with common sequences')
    parser.add_argument('-e','--extended', dest='phylip_type',
                        default=None, action='store_true',
                        help='If the phylip dataset is in the extended format, use this [individual names are still limited to a max of 15 characters].')
    parser.add_argument('-n','--num_loci', dest='num_loci',
                        default=1, action='store', type=int,
                        help='number of loci')
    parser.add_argument('-s','--species_name', dest='species_name',
                        default=None, action='store',type=str,
                        help='give a name to be used in nexus file as taxpartition name')
    parser.add_argument('-w','--write_myfile', dest='write_myfile',
                        action='store', default='myfile', type=str,
                        help='the name for the concatenated nexus-file/file of common sequences in all loci')
    parser.add_argument('-f','--folder', dest='folder',
                        action='store', default='loci', type=str,
                        help='the folder that contains the data (loci in separate text files called "loci0.txt", "loci1.txt", ...). The default is "loci"')
                        
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args
    
    
    
    
def read_data(current, folder, num_loci,prefix='locus',suffix='.txt'):
    labels=[]
    sequences=[]
    varsitess=[]
    for i in range(num_loci):
        locus_i = prefix+str(i)+suffix
        locus = os.path.join(current,folder,locus_i)
        label,sequence,varsites = phylip.readData(locus, ttype)
        if DEBUG:
            print(f"label = {label}")
            print(f"sequence = {sequence}")
        labels.append(label)
        sequences.append(sequence)
        varsitess.append(varsites)
    return [labels,sequences,varsitess]
    

    
    
def concatenate_sequences(myloci, prefix='locus',suffix='.txt'):
    labels, sequences, varsitess = myloci
    if DEBUG:
        for i in range(len(labels)):
            print(f"labels[{i}] = {labels[i]}")
    
    labels = [[x.lower() for x in label] for label in labels]
    print(f"\nlen(labels) = {[len(l) for l in labels]}")
    
    common = sorted(set.intersection(*map(set, labels)) , key = labels[0].index)   # common element extraction from N lists
    print(f"\ncommon = {common}")
    
    counting = [min([l.count(item) for l in labels]) for item in common]
    print(f"\ncounting = {counting}")
    
    numsites= [len(sequences[i][0]) for i in range(len(sequences))]
            
    dir = folder_copy
    if os.path.exists(dir):
        shutil.rmtree(dir)
    shutil.copytree(folder, dir )
    
    
    if max(counting)<=1:    # in the case in each locus we DO NOT have individuals with the same name ( like Scott Astralia bird dataset)
        common_indxs=[]
        for i in range(len(labels)):
            mylabel = labels[i]
            indxs = [mylabel.index(j) for j in common]
            common_indxs.append(indxs)
        if DEBUG:
            print(f"\ncommon_indxs = {common_indxs}")

        label=common
        common_sequences = [[sequences[i][j] for j in common_indxs[i]] for i in range(len(common_indxs))]
        
        
        if DEBUG:
            print(f"\ncommon_sequences = {common_sequences}")
            print(f"\nnumsites = {numsites}")
            
        for i in range(len(common_sequences)-1):
            sequence = [a+b for a,b in zip(common_sequences[i],common_sequences[i+1])]
            common_sequences[i+1]= sequence
        print(f"len(sequences) = {len(sequences)}")
        print(f"len(common_sequences) = {len(common_sequences)}")
        
        
        res_list = [[item for item in label if item not in common] for label in labels]
        if DEBUG:
            print(f"len(res_list) ={[len(item) for item in res_list]}")
        # in new folder, in each locus, remove sequences other than "common":
        for i in range(num_loci):
            locus_i = prefix+str(i)+suffix
            locus = os.path.join(current,folder_copy,locus_i)
            if DEBUG:
                print(f"i ={i} \nlocus = {locus}")
            with open(locus, 'r') as fp:
                lines = fp.readlines()
            
            i_numind,i_numsites, *rest = (lines[0]).split()
            if DEBUG:
                print(f"i_numind = {i_numind}, i_numsites = {i_numsites}")
            with open(locus, 'w') as fp:
                lines[0] = str(len(common)) +' '+ i_numsites+'\n'
                #check if any of the items of the list "common" is in a line of the file,if not, then remove that line
                for number, line in enumerate(lines):
                    if not any(word in line.lower() for word in res_list[i]):
                        fp.write(line)
            
    else:     # in the case in each locus we  have individuals with the same name ( like Fasta bird dataset)
        indxs = [[np.where(np.array(label) == c)[0].tolist() for c in common] for label in labels]
        
        remain_indxs=[]
        for indx in indxs:
            remain=[]
            for i in range(len(indx)):
                remain.extend(indx[i][:counting[i]])
            remain_indxs.append(remain)
        if DEBUG:
            print(f"\nremain_indxs = {remain_indxs}")
        common_num = sum(counting)
        if DEBUG:
            print(f"\ncommon_num = {common_num}")
        
        # in new folder, in each locus, remove sequences other than "common":
        for i in range(num_loci):
            locus_i = prefix+str(i)+suffix
            locus = os.path.join(current,folder_copy,locus_i)
            if DEBUG:
                print(f"i ={i} \nlocus = {locus}")
            with open(locus, 'r') as fp:
                lines = fp.readlines()
            i_numind,i_numsites, *rest = (lines[0]).split()

            if DEBUG:
                print(f"i_numind = {i_numind}, i_numsites = {i_numsites}")
            with open(locus, 'w') as fp:
                fp.write(str(common_num) +' '+ i_numsites+'\n')
                for number in range(1,len(lines)):
                    if number in remain_indxs[i]:
                        fp.write(lines[number])
        label= [labels[0][i] for i in remain_indxs[0]]
        if DEBUG:
            print(f"len(label) = {len(label)}")
            print(f"label = {label}")
        
        common_sequences=[]
        for i, loc_labels_num in enumerate(remain_indxs ):
            loc_seq = [sequences[i][j] for j in loc_labels_num]
            common_sequences.append(loc_seq)
        print(f"len(common_sequences) = {len(common_sequences)}")
        if DEBUG:
            print(f"common_sequences = {common_sequences}")
        
        for i in range(num_loci-1):
            sequence = [a+b for a,b in zip(common_sequences[i],common_sequences[i+1])]
            common_sequences[i+1]= sequence
        #sys.exit()
    return label, sequence, numsites


def write_files(label, sequence, numsites, letters, merging_nums, species_name, myfile):
    ntax= len(label)
    nchar = len(sequence[0])
    with  open(myfile+'.nex', "w") as f:
        f.write('#nexus\n\n')
        f.write('begin data;\n')
        f.write('         dimensions ntax={} nchar={};\n'.format(ntax, nchar))
        f.write('         format datatype=dna gap=-;\n')
        f.write('         matrix\n')
        for i in range(ntax):
            f.write(f'{label[i]:<20}{"".join(map(str, sequence[i]))}\n')
        f.write(';\nend;\nbegin sets;\n')
        f.write('    charpartition loci = locus0:1-{},\n'.format(numsites[0]))
        count = numsites[0]+1
        for i in range(1,len(numsites)-1):
            f.write('{}locus{}:{}-{},\n'.format(" "*25, i, count, count+numsites[i]-1))
            count +=numsites[i]
        i=len(numsites)-1
        f.write('{}locus{}:{}-{};\n'.format(" "*25, i, count, count+numsites[i]-1))
                    
        f.write('end;\n\n[ The taxpartition defined below assigns each of the snake samples to the appropriate species. ]\nbegin sets;\n')
        f.write('    taxpartition {} =\n'.format(species_name))
        count=1
        for i in range(len(letters)-1):
            myname=letters[i]
            f.write('{}{} : {}-{},\n'.format(" "*8,myname,count, count+len(merging_nums[i])-1))
            count +=len(merging_nums[i])
        i=len(letters)-1
        myname=letters[i]
        f.write('{}{} : {}-{};\n'.format(" "*8,myname,count, count+len(merging_nums[i])-1))
        f.write('end;\n')
    f.close()
    
    with  open(myfile, "w") as f:
        f.write('{} {}\n'.format(ntax, nchar))
        for i in range(ntax):
            f.write(f'{label[i]:<15}{"".join(map(str, sequence[i]))}\n')
    f.close()
    

if __name__ == '__main__':

    args = myparser() # parses the commandline arguments
    
    phylip_type = args.phylip_type
    if phylip_type:
        ttype = 'EXTENDED'
        filetype = 'RelPHYLIP'
    else:
        ttype = 'STANDARD'
        filetype = 'PHYLIP'
        
    num_loci = args.num_loci
    species_name = args.species_name
    myfile = args.write_myfile
    folder = args.folder
    folder_copy = folder+'_copy'

    current = os.getcwd()
    myloci = read_data(current, folder, num_loci, 'locus', '.txt')
    
    label, sequence, numsites = concatenate_sequences(myloci)

    if DEBUG:
        print(f"\nsequence = {sequence}")
    print(f"\nlabel = {label}")
    print(f"\nnumsites = {numsites}")

    

    letters = list(dict.fromkeys([x for x in label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"
    print(f"\nletters = {letters}")
        
    merging_nums=[]
    for item in letters:
        merg_indxs= [label.index(i) for i in label if item in i]
        merging_nums.append(merg_indxs)
        
    print(f"\nmerging_nums = {merging_nums}")
        
    write_files(label, sequence, numsites, letters, merging_nums, species_name, myfile)
    
                

    
