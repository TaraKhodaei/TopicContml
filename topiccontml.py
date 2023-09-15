#!/usr/bin/env python
#Import Libraries

# import numpy for matrix operation
import numpy as np

# Importing Gensim
import gensim
from gensim import corpora

# to generate dictionary of unique tokens
from gensim.corpora import Dictionary

#LDA model
from gensim.models import LdaModel, LdaMulticore

# to suppress warnings
from warnings import filterwarnings
filterwarnings('ignore')

#plotting tools
import matplotlib.pyplot as plt

#To extract my sequences from myfile.phy
import phylip

#To calculate distance of trees
import dendropy
from dendropy.calculate import treecompare


import os
import sys
import itertools

current = os.getcwd()


DEBUG = False
if DEBUG:
    #To print progress of the training procedure on screen
    import logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.NOTSET)
    
    
def myparser():
    import argparse
    parser = argparse.ArgumentParser(description='generates topic frequencies for the dataset to be applied as an input to CONTML to generate the phylogeny')
    parser.add_argument('-e','--extended', dest='phylip_type',
                        default=None, action='store_true',
                        help='If the phylip dataset is in the extended format, use this [individual names are still limited to a max of 15 characters].')
    parser.add_argument('-gt','--gaps_type', dest='gaps_type',
                        default=None, action='store',type=str,
                        help='String "rm_row": removes gaps(-) in each sequence (by row). String "rm_col": removes site columns thathave at least one gap(-). Without this option the gaps(-) are included.')
    parser.add_argument('-m','--merging', dest='merging',
                        default=None, action='store',type=int,
                        help='Merge sequences that start with the same 3 letters [e.g. population or species labels].')
    parser.add_argument('-kr','--kmer_range', dest='kmer_range',
                        default='2,10,2', action='store',
                        help='range of kmers extraction, lowerbound,max+1,step [for example: 2,10,2 leads to non overlapping k-mers: 2,4,6,8')
    parser.add_argument('-kt','--kmer_type', dest='kmers_type',
                        default='not_overlap', action='store',type=str,
                        help='default "not_overlap": extract kmers without overlapping. String "not_overlap": extract kmers with overlapping.')
    parser.add_argument('-f','--folder', dest='folder',
                        action='store', default='loci', type=str,
                        help='the folder that contains the data (loci in separate text files called "locus0.txt", "locus1.txt", ...). The default is "loci"')
    parser.add_argument('-n','--num_loci', dest='num_loci',
                        default=1, action='store', type=int,
                        help='number of loci')
    parser.add_argument('-sd','--siminfile_diverge_time', dest='sim_diverge',
                        default=None, action='store',
                        help='To do siminfile analysis for the folder with the given float number')
    parser.add_argument('-b','--bootstrap', dest='bootstrap',
                        default=0, action='store', type=int,
                        help='number of bootstrap replicates')
    parser.add_argument('-bt','--bootstrap_type', dest='bootstrap_type',
                        default='kmer', action='store',type=str,
                        help='default "kmer": do the bootsrap by randomly choosing  x kmers in each document of x kmers. String "seq": do the bootsrap by randomly choosing  x columns  of aligned sequences with the same length of x ("seq" works only in the ccase the sequences have the same lengths)')
    
                            
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args
    
    
#====================================================================================
#Tokenize Documents
def tokenize(seq_list, k, gaps_type, kmers_type):  #gaps:rm_row/rm_com     #kmers_type:overlap/not_overlap
    if gaps_type == 'rm_row':
        seq_list = [seq.replace('-', '') for seq in seq_list]
        
    elif gaps_type == 'rm_col':
        # String List to Column Character Matrix Using zip() + map()
        seq_matrix = np.transpose(list(map(list, zip(*seq_list))))
        #remove columns if the column has an element '-'
        seq_matrix_cleaned = np.array([col for col in seq_matrix.T if '-' not in col]).T
        #Convert List of lists to list of Strings again
        seq_list = [''.join(ele) for ele in seq_matrix_cleaned]
        
    if DEBUG:
        np.savetxt ('seq_list', seq_list,  fmt='%s')
        
    #docs: list of lists of each document kemers
    docs = []
    if kmers_type == 'overlap':
        for seq in seq_list:
            doc=[]
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                doc.append(kmer)
            docs.append(doc)
    elif kmers_type == 'not_overlap':
        for seq in seq_list:
            doc = [seq[i:i+k] for i in range(0, len(seq), k)]
            docs.append(doc)
    return docs


#====================================================================================
def read_data(current, folder, num_loci,prefix='locus',suffix='.txt'):
    labels=[]
    sequences=[]
    varsitess=[]
    for i in range(num_loci):
        #===================== Extract labels and sequences ====================
        locus_i = prefix+str(i)+suffix
        locus = os.path.join(current,folder,locus_i)        
        label,sequence,varsites = phylip.readData(locus, ttype)
        if DEBUG:
            print(f"len(label) = {len(label)}")
            print(f"len(sequence) = {len(sequence)}")
        labels.append(label)
        sequences.append(sequence)
        varsitess.append(varsites)
    return [labels,sequences,varsitess]



#====================================================================================
def topicmodeling(bootstrap, myloci, num_loci, num_topics, chunksize , passes , iterations , eval_every ):
    labels, sequences, varsitess = myloci
    topics_loci = []
    letters_loci = []
    for i in range(num_loci):
        print(f'\n~~~~~~~~~~~~~~~~~~~~~ locus {i} ~~~~~~~~~~~~~~~~~~~~~~~\n')
        label = labels[i]        
        sequence = sequences[i]
        varsites = varsitess[i]
        #========================= Extract k-mers ==============================
        docs=[]
        for k in range(kmer_range_list[0],kmer_range_list[1],kmer_range_list[2]):
            tokenize_k = tokenize(sequence, k, gaps_type, kmers_type)
            if len(docs)>0:
                docs =[docs[i]+tokenize_k[i] for i in range(len(tokenize_k))]
            else:
                docs = tokenize_k
        if bootstrap == 'kmer' and nbootstrap>0:
            #sys.exit()
            print(f"bootstrapping k-mers")
            for i in range(len(docs)):
                docs[i] = np.random.choice(docs[i],size=len(docs[i]),replace=True)
            
        #count number of all words in docs
        count = 0
        for doc in docs:
            count += len(doc)
        print(f"Number of all words in docs = {count}")
        
        
        letters = label
        if DEBUG:
            print(f"len(letters) = {len(letters)} ")
            print(f"letters = {letters} ")

        
        #====================== k-mers after merging ===========================

        if merging:

            letters = list(dict.fromkeys([x[:merging] for x in label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"
            if DEBUG:
                print(f"letters = {letters} ")
        
            merging_num=[]
            for item in letters:
                merg_indxs= [label.index(i) for i in label if item in i]
                merging_num.append(merg_indxs)
                
            merging_nums=[len(i) for i in merging_num]
            print(f"merging_nums = {merging_nums}")
            
            docs_merged=[]
            j=0
            for num in merging_nums:
                doc_merged = list(itertools.chain.from_iterable(docs[j:j+num]))
                docs_merged.append(doc_merged)
                j +=num
            print(f"len(docs_merged) = {len(docs_merged)} ")
            
            
            #count number of all words in docs
            count = 0
            for doc in docs_merged:
                count += len(doc)
                
            if DEBUG:
                print(f"Number of all words in docs_merged = {count}")
            
            docs = docs_merged
#        sys.exit()    #debug
        #================= Dictionary of Unique Tokens ===========================
        
        dictionary = Dictionary(docs)
        print(f"\nDictionary:\n{dictionary}")
        
        #=============== Filtering: remove rare and common tokens ================
        dictionary.filter_extremes(no_below=2, no_above=0.5)
        print(f"\nDictionary after filtering:\n{dictionary}")
        
        #================ Vectorize data: Bag-of-words ===========================
        corpus = [dictionary.doc2bow(doc) for doc in docs]
        print(f'\nNumber of documents: {len(corpus)}')
        print(f'Number of unique tokens: {len(dictionary)}')
        
        #===================== LDA Model: Training ===============================

        # Make a index to word dictionary.
        temp = dictionary[0]
        id2word = dictionary.id2token

        model = LdaModel(corpus=corpus, id2word=id2word, chunksize=chunksize, \
                           alpha=1, eta='auto', \
                           iterations=iterations, num_topics=num_topics, \
                           passes=passes, eval_every=eval_every, minimum_probability=0, update_every=5 )

        #================== Print topics/words frequencies ======================
        if DEBUG:
            for idx, topic in model.print_topics(-1):
                print("Topic: {} \nWords: {}\n".format(idx, topic))
                
        #============== Assigning the topics to the documents ===================
        docs_tuples = model[corpus]
                
        #=================== topics list for current locus ======================
        topics = []
        for num, doc in enumerate(docs_tuples):
            first=[]
            second=[]
            for tuple_i in doc:
                first.append(tuple_i[0])
                second.append(tuple_i[1])
                
            topics_freq=[]
            for i in range(num_topics):
                topics_freq.append(second[first.index(i)])

            topics.append(topics_freq)

        topics_loci.append(topics)
        letters_loci.append(letters)
        
        
        common_letters = list(set.intersection(*map(set, letters_loci)))
        if DEBUG:
            print(f"common_letters = {common_letters}")
        common_letters_loci_indx = [[l.index(c) for c in common_letters] for l in letters_loci]
        if DEBUG:
            print(f"common_letters_loci_indx = {common_letters_loci_indx}")
        common_topics_loci = [[l[i] for i in common_letters_loci_indx[num]] for num,l in enumerate(topics_loci)]
        if DEBUG:
            print(f"common_topics_loci = {common_topics_loci}")
        

    return common_letters, common_topics_loci

    

#====================================================================================
def infile(topics_loci, letters, num_loci):
    print(f"num_loci3 = {num_loci}")
    letters_limited = [x[:10] for x in letters]      #CONTML accepts population names lass than 10 letters
    num_pop= len(topics_loci)
    with  open("infile", "w") as f:
        f.write('     {}    {}\n'.format(num_pop, num_loci))
        f.write('{} '.format(num_topics)*num_loci)
        f.write('\n')
        for i in range(num_pop):
            myname = letters[i]
            f.write(f'{myname:<15}{" ".join(map(str, topics_loci[i]))}\n')
    f.close()

#====================================================================================
def run_contml(infile):
    os.system('rm outfile outtree')
    contmljumble = RANDOMSEED
    contmltimes  = 10
    with open('contmlinput','w') as f:
        contmlinput = f'g\nj\n{contmljumble}\n{contmltimes}\ny'
        #contmlinput = f'g\n'
        f.write(contmlinput)
    os.system(f'cat contmlinput | {PROGRAMPATH}contml2')
#    os.system(f'cat contmlinput | {PROGRAMPATH}contml')
    
    #read the outtree file
    with open('outtree', 'r') as f:
        tree = f.read().replace('\n', ' ')
    return tree

#====================================================================================
def simulation(current, folder):
    diverge_time = float(sim_diverge)       #diverge_time=0.0/0.01/0.05/0.1/0.2
    tree_newick = '((Arb:1,Fra:1):1,Kre:1,Rom:1);'    #We need just topology to compare
    tns = dendropy.TaxonNamespace()
    true_tree = dendropy.Tree.get(data=tree_newick,schema="newick",taxon_namespace=tns)
    true_tree.encode_bipartitions()
        
    simfiles_folder = os.path.join(current,folder)
    print(simfiles_folder)
    current = simfiles_folder
    siminfiles_trees=[]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    num_iter =[2,5,10,20,50,100]
    sims_num = list(range(0,100))
    count_equaltrue=[]
    for num in num_iter:
        num_loci = num
        distances =[]
        for sim_num in sims_num:
            siminfile_sl ='siminfile_'+str(sim_num)+'_100_'+str(diverge_time)+'_100'
            folder = os.path.join(simfiles_folder,siminfile_sl)
            #sys.exit()
            loci = read_data(current, folder, num_loci, 'locus', '.txt')
            letters, topics_loci = topicmodeling('NoBootstrap',loci, num_loci, num_topics, chunksize , passes , iterations , eval_every )
            
            #for each locus remove last column from topic matrix
            topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
            
            
            #concatenation of the topics for all loci
            topics_loci_concatenated = topics_loci_missingLast[0]
            for i in range(1,len(topics_loci_missingLast)):
                topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]

            #generate infile
            infile(topics_loci_concatenated, letters, num_loci)
            
            #run CONTML
            ourtree = run_contml(infile)
            print(ourtree)
            
            our_tree = dendropy.Tree.get(data=ourtree,schema="newick",taxon_namespace=tns)
            our_tree.encode_bipartitions()
            distance=treecompare.unweighted_robinson_foulds_distance(true_tree, our_tree, is_bipartitions_updated=True)
            distances.append(distance)
            #Figtree
            #os.system(f"{PROGRAMPATH}figtree outtree")
        count_equaltrue.append(distances.count(0))
            
    np.savetxt('trueAgreement_simulation{}'.format(diverge_time), count_equaltrue,  fmt='%s')
    print(f"\nResults of agreement with true tree using RF-distanc is witten to the file 'trueAgreement_simulation{diverge_time}'")



#====================================================================================
# a single analysis that shows the tree using figtree with show=True
def single_run(show=True):
    loci = read_data(current, folder, num_loci, 'locus', '.txt')
    letters, topics_loci = topicmodeling('NoBootstrap',loci, num_loci, num_topics, chunksize , passes , iterations , eval_every )
    
    #for each locus remove last column from topic matrix
    topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
    
    
    #concatenation of the topics for all loci
    topics_loci_concatenated = topics_loci_missingLast[0]
    for i in range(1,len(topics_loci_missingLast)):
        topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
        #print(f'topics_loci_concatenated =\n{topics_loci_concatenated}')
        
    #generate infile
    infile(topics_loci_concatenated, letters, num_loci)
    
    #run CONTML
    ourtree = run_contml(infile)
    print(ourtree)
    if show:
        #Figtree
        os.system(f"{PROGRAMPATH}figtree outtree")
    return ourtree

        
#====================================================================================
def bootstrap_sequences(sequences):
    newsequences = []
    for locus in sequences:
        nind = len(locus)
        nsites = len(locus[0])
        pick = np.random.randint(0,nsites,nsites)
        locusnewseq=[]
        for ni in range(nind):
            newseq = "".join([locus[ni][i] for i in  pick])
            locusnewseq.append(newseq)
        newsequences.append(locusnewseq)
    return newsequences
    
#====================================================================================
def bootstrap_run(bootstrap):
    count_boot = 1
    outtrees=[]
    loci = read_data(current, folder, num_loci, 'locus', '.txt')
    for bi in range(nbootstrap):
        print(f"TEST===> bootstrap ={bootstrap}")
        labels,sequences,varsitess = loci
        if not bootstrap=='kmer':      #in case of "seq"
            bsequences = bootstrap_sequences(sequences)
            bloci = [labels,bsequences,varsitess]
        else:
            bloci = loci
        count_boot += 1
        #sys.exit()      
        letters, topics_loci = topicmodeling(bootstrap,bloci, num_loci, num_topics, chunksize , passes , iterations , eval_every )
        #for each locus remove last column from topic matrix
        topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
        #concatenation of the topics for all loci
        topics_loci_concatenated = topics_loci_missingLast[0]
        for i in range(1,len(topics_loci_missingLast)):
            topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
        #generate infile
        infile(topics_loci_concatenated, letters, num_loci)
        #run CONTML
        ourtree = run_contml(infile)
        outtrees.append(ourtree)
    return outtrees



#====================================================================================

if __name__ == "__main__":
    args = myparser() # parses the commandline arguments
    gaps_type = args.gaps_type
    merging = args.merging
    sim_diverge = args.sim_diverge
    kmers_type = args.kmers_type
    num_loci = args.num_loci
    folder = args.folder
    bootstrap_type = args.bootstrap_type

    nbootstrap = args.bootstrap
    if nbootstrap == 0:
        bootstrap = 'none'
    else:
        bootstrap = bootstrap_type        #bootstrap = 'kmer' or 'seq'
    
    phylip_type = args.phylip_type
    if phylip_type:
        ttype = 'EXTENDED'
        filetype = 'RelPHYLIP'
    else:
        ttype = 'STANDARD'
        filetype = 'PHYLIP'

    kmer_range_list = list(map(int, args.kmer_range.split(',')))
    
    num_topics= 5
    chunksize = 20
    passes = 50
    iterations = 1000
    eval_every = 1


    PROGRAMPATH = '/Users/tara/bin/'
    # generates a random number seed for jumble in contml
    RANDOMSEED  = np.random.randint(1,2**16,size=1)[0]
    if RANDOMSEED % 2 == 0:
        RANDOMSEED += 1
     

    #=================== main: if siminfile analysis ======================
    if sim_diverge is not None:
        simulation(current, folder)
    #===================  main: if NOT siminfile analysis  ======================
    elif bootstrap!= 'none':
        ourtree = single_run(False)
        with open('best.tre','w') as btrees:
            btrees.write(ourtree+'\n')
            
        outtrees = bootstrap_run(bootstrap)
        with open('bootstrap_replicates.tre','w') as btrees:
            for tr in outtrees:
                btrees.write(tr+'\n')
        os.system(f"sumtrees.py --decimals=0 --percentages --output-tree-filepath=bootstrap_target_best.tre --target=best.tre bootstrap_replicates.tre")
        os.system(f"sumtrees.py --decimals=0 --percentages --output-tree-filepath=bootstrap.tre bootstrap_replicates.tre")
    else:
        single_run()
