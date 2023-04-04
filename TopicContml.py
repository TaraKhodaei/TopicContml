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


import os
import itertools

current = os.getcwd()
#print(current)

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
                        help='If the phylip dataset is in the extended format, use this.')
    parser.add_argument('-gt','--gaps_type', dest='gaps_type',
                        default=None, action='store',type=str,
                        help='String "rm_row": removes gaps(-) in each sequence by the row. String "rm_col": romoves the column if there is at least one gap(-) in that column. Otherwise, it keeps the gaps(-) ')
    parser.add_argument('-m','--merging', dest='merging',
                        default=None, action='store_true',
                        help='It merges sequences with the same population.')
    parser.add_argument('-kr','--kmer_range', dest='kmer_range',
                        default='2,10,2', action='store',
                        help='range of kmers extraction')
    parser.add_argument('-kt','--kmer_type', dest='kmers_type',
                        default='not_overlap', action='store',type=str,
                        help='default "not_overlap": extract kmers without overlapping. String "not_overlap": extract kmers with overlapping.')
    parser.add_argument('-n','--num_loci', dest='num_loci',
                        default=1, action='store', type=int,
                        help='number of loci')

                            

    args = parser.parse_args()
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
def topicmodeling(num_loci, num_topics, chunksize , passes , iterations , eval_every ):
    topics_loci=[]
    for i in range(num_loci):
        print(f'\n~~~~~~~~~~~~~~~~~~~~~ locus {i} ~~~~~~~~~~~~~~~~~~~~~~~\n')
                
        #===================== Extract labels and sequences ====================
        locus_i = 'locus'+str(i)+'.txt'
        locus = os.path.join(current,'loci',locus_i)
        
        label,sequence,varsites = phylip.readData(locus, ttype)
        #print(f"label = {label}")
        #print(f"sequence = {sequence}")
        
        #========================= Extract k-mers ==============================
        docs=[]
        for k in range(kmer_range_list[0],kmer_range_list[1],kmer_range_list[2]):
            tokenize_k = tokenize(sequence, k, gaps_type, kmers_type)
            if len(docs)>0:
                docs =[docs[i]+tokenize_k[i] for i in range(len(tokenize_k))]
            else:
                docs = tokenize_k
                
        #count number of all words in docs
        count = 0
        for doc in docs:
            count += len(doc)
        print(f"Number of all words in docs = {count}")
        
        
        letters = label
        
        #====================== k-mers after merging ===========================
        if merging:
            #extraxt distinct sequence names based on first 3 letters
            letters = list(dict.fromkeys([x[:3] for x in label ]))  #remove duplicates from a list, while preserving order using"dict.fromkeys" insted of "list(set)"
            print(letters)
            
            merging_nums=[]
            for item in letters:
                merg_indxs= [label.index(i) for i in label if item in i]
                merging_nums.append(merg_indxs)
                
            merging_nums=[len(i) for i in merging_nums]
            print(f"merging_nums = {merging_nums}")
            
            docs_merged=[]
            j=0
            for num in merging_nums:
                doc_merged = list(itertools.chain.from_iterable(docs[j:j+num]))
                docs_merged.append(doc_merged)
                j +=num

            #count number of all words in docs
            count = 0
            for doc in docs_merged:
                count += len(doc)
            print(f"Test ===> Number of all words in docs_merged = {count}")
            
            docs = docs_merged
            
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
        delta_type='dirichlet'
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
            
        #print(f'\nTopics :\n{topics}')
        topics_loci.append(topics)
    return letters, topics_loci

#====================================================================================
def infile(topics_loci, letters):
    letters_limited = [x[:10] for x in letters]      #CONTML accepts population names lass than 10 letters
    num_pop= len(topics_loci)
    with  open("infile", "w") as f:
        f.write('     {}    {}\n'.format(num_pop, num_loci))
        f.write('{} '.format(num_topics)*num_loci)
        f.write('\n')
        for i in range(num_pop):
            myname = letters[i]
            f.write(f'{myname}{" "*(11-len(myname))}{" ".join(map(str, topics_loci[i]))}\n')
    f.close()

#====================================================================================
def run_contml(infile):
    os.system('rm outfile outtree')
    os.system('echo "y" | /Users/tara/bin/contml')
    
    #read the outtree file
    with open('outtree', 'r') as f:
        tree = f.read().replace('\n', ' ')
    return tree



#====================================================================================

if __name__ == "__main__":
    
    args = myparser() # parses the commandline arguments
    
    gaps_type = args.gaps_type
    
    merging = args.merging
    
    kmers_type = args.kmers_type
    
    num_loci = args.num_loci
    
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
    
    letters, topics_loci = topicmodeling(num_loci, num_topics, chunksize , passes , iterations , eval_every )

    #for each locus remove last column from topic matrix
    topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]


    #concatenation of the topics for all loci
    topics_loci_concatenated = topics_loci_missingLast[0]
    for i in range(1,len(topics_loci_missingLast)):
        topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
    #print(f'topics_loci_concatenated =\n{topics_loci_concatenated}')



    #generate infile
    infile(topics_loci_concatenated, letters)


    #run CONTML
    ourtree = run_contml(infile)
    print(ourtree)


    #Figtree
    os.system("/Users/tara/bin/figtree outtree")
