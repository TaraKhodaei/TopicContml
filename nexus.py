#!/usr/bin/env python
# coding: utf-8

import phylip
import os
import sys
import shutil
import numpy as np
from itertools import chain

DEBUG=False


    
def myparser():
    import argparse
    parser = argparse.ArgumentParser(description='generates nexus file of concatenated loci with common sequences')
    parser.add_argument('-e','--extended', dest='phylip_type',
                        default=None, action='store_true',
                        help='If the phylip dataset is in the extended format, use this [individual names are still limited to a max of 15 characters].')
    parser.add_argument('-nl','--num_loci', dest='num_loci',
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
    parser.add_argument('-m','--merging', dest='merging',
                        default=None, action='store',type=int,
                        help='Merge sequences that start with the same 3 letters [e.g. population or species labels].')
    parser.add_argument('-t','--type', dest='nexus_type',
                        default='total', action='store', type=str,
                        help='how to generate the nexus file. If "common", then it generates a merged nexus file of sequences with common labels in all loci. Default "total" generates a merged nexus file of sequence with totol labels in all loci such that in each locus it adds sequences with all gaps if thaat sequence is not in that locus.')
                        
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args
    
    
    
#-------------------------------------------------------
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
    

    
#-------------------------------------------------------
#def concatenate_sequences(myloci, prefix='locus',suffix='.txt'):
def concatenate_sequences(myloci, nexus_type, prefix='locus',suffix='.txt'):
    labels, sequences, varsitess = myloci
    

    if DEBUG:
        for i in range(len(labels)):
            print(f"labels[{i}] = {labels[i]}")
    

    labels = [[x.lower() for x in label] for label in labels]
    print(f"\nlen(labels) = {[len(l) for l in labels]}")

    
    common = sorted(set.intersection(*map(set, labels)) , key = labels[0].index)   # common element extraction from N lists
    print(f"\n\ncommon = {common}")
    print(f"len(common) = {len(common)}")
    
    
    counting = [min([l.count(item) for l in labels]) for item in common]
    print(f"\ncounting = {counting}")
    
    numsites= [len(sequences[i][0]) for i in range(len(sequences))]
            
    dir = folder_copy
    if os.path.exists(dir):
        shutil.rmtree(dir)
    shutil.copytree(folder, dir )
    
    #+++++++++++++++++++++++++++  loci DO NOT have individuals with the same name  +++++++++++++++++++++++++++
    
    if max(counting)<=1:    # in the case in each locus we DO NOT have individuals with the same name ( like Scott Astralia bird dataset)
        print(f"\n\nmax(counting)=0 ===> loci DO NOT have individuals with the same name")
        common_indxs=[]
        for i in range(len(labels)):   #length of num_loci
            locuslabels = labels[i]
#            print(f"\nTEST:\nlen(locuslabels) = {len(locuslabels)}")
            indxs = [locuslabels.index(j) for j in common]
#            print(f"TEST:\nlen(indxs) = {len(indxs)}")
            common_indxs.append(indxs)
#        if DEBUG:
        print(f"\ncommon_indxs = {common_indxs}")

        label=common
#        common_sequences = [[sequences[i][j] for j in common_indxs[i]] for i in range(len(common_indxs))]
        common_sequences = [[sequences[i][j] for j in common_indxs[i]] for i in range(num_loci)]
        
#        for i in range(num_loci):   #test  OK
#            print(f"\n{common_sequences[i][0]}")
        
        print(f"\n\nTEST:\nlen(common_sequences) = {len(common_sequences)}")
        
        if DEBUG:
            print(f"\ncommon_sequences = {common_sequences}")
            print(f"\n\nnumsites = {numsites}")
        
        #-----------------
#        print(f"\n\nlabels = {labels}")
        flat_labels = list(chain.from_iterable(labels))
        total_labels = list(set(','.join(flat_labels).split(',')))        # total labels in all loci
        #if DEBUG:
#        print(f"\n\nflat_labels = {flat_labels}")
        print(f"\n\ntotal_labels = {total_labels}")
        print(f"len(total_labels) = {len(total_labels)}")
        
        uncommon = sorted(list(set(total_labels) - set(common) ))
        #if DEBUG:
        print(f"\n\nuncommon = {uncommon}")
        print(f"len(uncommon) = {len(uncommon)}")
        
        uncommon_sequences = []
        for i in range(num_loci):
            seq_list=[]
            for labl in uncommon:
                if labl in labels[i]:
                    l = labels[i].index(labl)      #Note: if it is unique
                    seq_list.append(sequences[i][l])
                else:
                    seq_fill ='?'*len(sequences[i][0])
                    seq_list.append(seq_fill)
            uncommon_sequences.append(seq_list)
            
#        for i in range(num_loci):   #test OK
#            print(f"\n{uncommon_sequences[i][0]}")
            
        if DEBUG:
    #        print(f"TEST:\n\nuncommon sequences for first locus = {uncommon_sequences[0]}")
            print(f"\n\n============== TEST================== ")
            print(f"len(uncommon sequences for first locus) = {len(uncommon_sequences[0])}")
            print(f"\nlen(uncommon sequences for second locus) = {len(uncommon_sequences[1])}")
            print(f"\nlen(uncommon sequences for third locus) = {len(uncommon_sequences[2])}")
            print(f"\n============== TEST================== ")
        
        #===========================================
        if nexus_type == "common":
            print(f"TEST : nexus_type ={nexus_type}")
            
            label=common
                    
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
        #===========================================
        elif nexus_type == "total":
            print(f"TEST : nexus_type ={nexus_type}")
            #creating merged sequences to do bootstrap using SVDquartets:
            
#            label=total        #First bug is here: it should be in order of common+uncommon not mixed labels !!!!!!
            label= common + uncommon

            
         
            total_seqs= [a + b for a, b in zip(common_sequences, uncommon_sequences)]   #common_sequences and uncommon_sequences are both length of num_loci--> in each locus: first common seqs + then uncommon seqs
            
            if DEBUG:
                print(f"\n\n============== TEST================== ")
                print(f"len(common_sequences) ={len(common_sequences)}")
                print(f"len(uncommon_sequences) ={len(uncommon_sequences)}")
                print(f"len(total_seqs) ={len(total_seqs)}")
                print(f"len(total_seqs[0]) for locus0 ={len(total_seqs[0])}")
                print(f"\n============== TEST================== ")
            
#            for i in range(len(total_seqs)-1):   #length of num_loci-1
            for i in range(num_loci-1):   #length of num_loci-1
                sequence = [a+b for a,b in zip(total_seqs[i],total_seqs[i+1])]
                total_seqs[i+1]= sequence


            #creating folder copy to do bootstrap using TopicContml:
            dir = folder_copy
            if os.path.exists(dir):
                shutil.rmtree(dir)
            shutil.copytree(folder, dir )
                    
            
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
                    
                res_labels = sorted(list(set(total_labels) - set(labels[0]) ))
#                seq_qap ='-'*int(i_numsites)
                seq_qap ='?'*int(i_numsites)
                
                    
                with open(locus, 'a') as fp:    #Open the file for 'append' rather than 'write'
                    #check if any of the items of the list "common" is in a line of the file,if not, then remove that line
                    #NOTE: labels here are all lowecase
                    for i in range(len(res_labels)):
                        fp.write(f'{res_labels[i]:<15}{"".join(map(str, seq_qap))}\n')
                        
                with open(locus, 'r') as fp:
                    data = fp.readlines()
                data[0] = str(len(total_labels)) +' '+ i_numsites+'\n'
                if DEBUG:
                    print(f"\n\ndata[1] ={data[1]} ")
                    print(f"data[1][0] ={data[1][0]} ")
                    print(f"data[1][1] ={data[1][1]} ")
                with open(locus, 'w' ) as fp:
                    fp.writelines(data)
                    
                    
                                        
                    
    #+++++++++++++++++++++++++++  loci DO have individuals with the same name  +++++++++++++++++++++++++++
    
    else:     # in the case in each locus we  have individuals with the same name ( like Fasta bird dataset)
        print(f"\n\nmax(counting)>0 ===> loci DO have individuals with the same name")
        indxs = [[np.where(np.array(label) == c)[0].tolist() for c in common] for label in labels]
        print(f"\nindxs = {indxs}")
                
        remain_indxs=[]
        for indx in indxs:
            remain=[]
            for i in range(len(indx)):
                remain.extend(indx[i][:counting[i]])
            remain_indxs.append(remain)
#        if DEBUG:
        print(f"\nremain_indxs = {remain_indxs}")
        common_num = sum(counting)
#        if DEBUG:
        print(f"\ncommon_num = {common_num}")
        
        #----------------------
        # in new folder, in each locus, remove sequences other than "common"
        for i in range(num_loci):
            locus_i = prefix+str(i)+suffix
            locus = os.path.join(current,folder_copy,locus_i)
#            if DEBUG:
#            print(f"i ={i} \nlocus = {locus}")
            with open(locus, 'r') as fp:
                lines = fp.readlines()
            i_numind,i_numsites, *rest = (lines[0]).split()

#            if DEBUG:
#            print(f"i_numind = {i_numind}, i_numsites = {i_numsites}")
            with open(locus, 'w') as fp:
                fp.write(str(common_num) +' '+ i_numsites+'\n')
                for number in range(1,len(lines)):
                    if number in remain_indxs[i]:
                        fp.write(lines[number])
        #----------------------
        
        label= [labels[0][i] for i in remain_indxs[0]]
#        if DEBUG:
        print(f"len(label) = {len(label)}")
        print(f"label = {label}")
        
        common_sequences=[]
        for i, loc_labels_num in enumerate(remain_indxs ):
            loc_seq = [sequences[i][j] for j in loc_labels_num]
            common_sequences.append(loc_seq)
        print(f"\nlen(common_sequences) = {len(common_sequences)}")
        print(f"len(common_sequences[0]) = {len(common_sequences[0])}")
#        if DEBUG:
#        print(f"common_sequences = {common_sequences}")
        
        for i in range(num_loci-1):
            sequence = [a+b for a,b in zip(common_sequences[i],common_sequences[i+1])]
            common_sequences[i+1]= sequence
        #sys.exit()
    return label, sequence, numsites






#-------------------------------------------------------
def write_files(label, sequence, numsites, letters, merging_nums, species_name, myfile):
    #Second bug is here: inmyfile.nex last paragraph note thaat "taxpartition birdspecies =" shuld be in the same order as appeard in file. Then it is better to write the seqs and labels with sorted one !!!!!!
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
    
    
    
#############################################################################
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
    nexus_type = args.nexus_type
    merging = args.merging

    current = os.getcwd()
    myloci = read_data(current, folder, num_loci, 'locus', '.txt')
    
#    label, sequence, numsites = concatenate_sequences(myloci)
    label, sequence, numsites = concatenate_sequences(myloci, nexus_type)

    if DEBUG:
        print(f"\nsequence = {sequence}")
    print(f"\nlabel = {label}")
    print(f"len(label) = {len(label)}")
    
    print(f"\n\nnumsites = {numsites}")
    
    
    
    
    
    #debuging:
    sort_label = sorted(label)
    print(f"\nsort_label = {sort_label}")
    sort_index = [i for i, x in sorted(enumerate(label), key=lambda x: x[1])]
    print(f"\nsort_index = {sort_index}")
    sort_sequence = [sequence[i] for i in sort_index]
#    print(f"\nsort_sequence = {sort_sequence}")

    
    
    
    
    

#    letters = list(dict.fromkeys([x for x in label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"

#    letters = list(dict.fromkeys([x[:merging] for x in label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"
    letters = list(dict.fromkeys([x[:merging] for x in sort_label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"      #debug
    print(f"\n\nletters = {letters}")
    print(f"len(letters) = {len(letters)}")
        
#    merging_nums=[]
#    for item in letters:
#        merg_indxs= [label.index(i) for i in label if item in i]
#        merging_nums.append(merg_indxs)
    merging_nums=[]    #debug
    for item in letters:
        merg_indxs= [sort_label.index(i) for i in sort_label if item in i]
        merging_nums.append(merg_indxs)
        
    print(f"\n\nmerging_nums = {merging_nums}")
    print(f"len(merging_nums) = {len(merging_nums)}")
        
#    write_files(label, sequence, numsites, letters, merging_nums, species_name, myfile)
    write_files(sort_label, sort_sequence, numsites, letters, merging_nums, species_name, myfile)      #debug
    
                

    
