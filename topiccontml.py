#!/usr/bin/env python
# MIT License
# (c) Tara Khodaei and Peter Beerli 2023
# github: khodaei
#Import Libraries
import numpy as np
import multiprocessing
from  tqdm import tqdm
from multiprocessing import Pool
import gensim
from gensim import corpora
from gensim.corpora import Dictionary
from gensim.models import LdaModel, LdaMulticore
from warnings import filterwarnings
filterwarnings('ignore')
import matplotlib.pyplot as plt
import phylip
import dendropy
from dendropy.calculate import treecompare
import os
import sys
import itertools
import time
from itertools import chain
from collections import Counter

#---------------------------------------------------
# most important global variables are defined in __main__
DEBUG = False
#---------------------------------------------------



def citations(options):
    print(" --------------------------------------------------------------------------------- ")
    print("| If you use this software for publications please cite these:                    |")
    print("|                                                                                 |")
    print("| Khodaei, M., Edwards, S. Beerli, P. (2023). Multilocus Phylogeny Estimation     |")
    print("|     Using Probabilistic Topic Modeling, Biorxiv doi: xxxx                       |")
    print("| Blei, D. M., Ng A. Y, and Jordan, M. I. (2003). Latent Dirichlet allocation.    |")
    print("|     Journal of machine Learning research, 3:993--1022                            |")
    print("| Felsenstein, J. (2005). PHYLIP (Phylogeny Inference Package) version 3.6.       |")
    print("|     Distributed by the author. Department of Genome Sciences, University        |")
    print("|     of Washington, Seattle. (https://phylipweb.github.io/phylip/)               |")
    print("| Řehůřek, R., and Sojka, P. (2010). Software framework for topic modelling with  |")
    print("|     large corpora. In proceedings of LREC 2010 Workshop on New Challenges       |")
    print("|     for NLP Frameworks, Valletta, Malta, pp.45--50.                              |")
    print("|     (http://is.muni.cz/publication/884893/en)                                   |")
    if bootstrap!= 'none':
        print("| Sukumaran, J. and Holder, M. T. (2010). DendroPy: a Python library for          |")
        print("|     phylogenetic computing, Bioinformatics, 26:1569--1571.                       |")
        print("|     (https://dendropy.org/)                                                     |")
    print(" --------------------------------------------------------------------------------- ")

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
                        default=None, action='store',  type=str,
                        help='the folder that contains loci data in separate text files called "locus0.txt", "locus1.txt", ...')
    parser.add_argument('-nf','--nexus_file', dest='nexus_file',
                        default=None, action='store',  type=str,
                        help='the nexus file that contains the multiloci data')
    parser.add_argument('-nl','--num_loci', dest='num_loci',
                        default=1, action='store', type=int,
                        help='number of loci')
    parser.add_argument('-sd','--siminfile_diverge_time', dest='sim_diverge',
                        default=None, action='store',
                        help='To do siminfile analysis for the folder with the given float number')
    parser.add_argument('-nb','--num_bootstrap', dest='num_bootstrap',
                        default=0, action='store', type=int,
                        help='number of bootstrap replicates')
    parser.add_argument('-bt','--bootstrap_type', dest='bootstrap_type',
                        default='kmer', action='store',type=str,
                        help='default "kmer": do the bootsrap by randomly choosing  x kmers in each document of x kmers. String "seq": do the bootsrap by randomly choosing  x columns  of aligned sequences with the same length of x ("seq" works only in the ccase the sequences have the same lengths)')
    parser.add_argument('-incl','--include', dest='include_file',
                        default=None, action='store',  type=str,
                        help='the include file contains a list of names that must be analyzed')
    parser.add_argument('-excl','--exclude', dest='exclude_file',
                        default=None, action='store',  type=str,
                        help='the exclude file contains a list of names that should not be analyzed')
    parser.add_argument('-force','--force', dest='force',
                        default=False, action='store_true',
                        help='this forces to use all species using an uninformative topicfrequency for missings')

    parser.add_argument('-show','--showtree', dest='showtree',
                        default=False, action='store_true',
                        help='uses figtree to show the tree')
    

    #gensim LDA arguments:
    parser.add_argument('-nt','--num_topics', dest='num_topics',
                        default=5, action='store', type=int,
                        help='Number of requested latent topics to be extracted from the training corpus. Defult value is 5 topics.')
    parser.add_argument('-i','--iterations', dest='iterations',
                        default=1000, action='store', type=int,
                        help='Maximum number of iterations through the corpus when inferring the topic distribution of a corpus. Defult value is 1000 iterations.')
    parser.add_argument('-p','--passes', dest='passes',
                        default=50, action='store', type=int,
                        help='Number of passes through the corpus during training. Defult value is 50.')
    parser.add_argument('-cs','--chunksize', dest='chunksize',
                        default=20, action='store', type=int,
                        help='Number of documents to be used in each training chunk. Defult value is 20.')
    parser.add_argument('-ee','--eval_every', dest='eval_every',
                        default=1, action='store', type=int,
                        help='Log perplexity is estimated every that many updates. Defult value is 1.')
    parser.add_argument('-ue','--update_every', dest='update_every',
                        default=5, action='store', type=int,
                        help='Number of documents to be iterated through for each update. Defult value is 5.')
    parser.add_argument('-al','--alpha', dest='alpha',
                        default='1',
                        help='a priori belief on document-topic distribution. It can be: (1) scalar for a symmetric prior over document-topic distribution, (2) 1D array of length equal to num_topics to denote an asymmetric user defined prior for each topic. (3) Alternatively default prior strings:"symmetric": a fixed symmetric prior of 1.0 / num_topics,"asymmetric": a fixed normalized asymmetric prior of 1.0 / (topic_index + sqrt(num_topics)),"auto":Learns an asymmetric prior from the corpus')
    parser.add_argument('-et','--eta', dest='eta',
                        default='auto',
                        help='a priori belief on topic-word distribution. It can be: (1) scalar for a symmetric prior over  topic-word distribution, (2) 1D array of length equal to num_words to denote an asymmetric user defined prior for each word, (3) matrix of shape (num_topics, num_words) to assign a probability for each word-topic combination. (4) Alternatively default prior strings:"symmetric": a fixed symmetric prior of 1.0 / num_topics,"auto": Learns an asymmetric prior from the corpus.')
                            
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args
    
    
#====================================================================================
#Tokenize Documents
def tokenize(seq_list, k, gaps_type, kmers_type):  #gaps:rm_row/rm_com     #kmers_type:overlap/not_overlap
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
def read_data(current, folder, num_loci,prefix='locus',suffix='.txt', include_names=set([]), exclude_names=set([])):
    labels=[]
    sequences=[]
    varsitess=[]
    seti = set(include_names)
    sete = set(exclude_names)
    real_loci = 0
    for i in range(num_loci):
        #Extract labels and sequences:
        locus_i = prefix+str(i)+suffix
        locus = os.path.join(current,folder,locus_i)        
        label,sequence,varsites = phylip.readData(locus, ttype)
        setl = set(label)
        iok = (setl.intersection(seti) == seti)
        if iok:
            real_loci += 1
            setle = setl.intersection(sete)
            setl  = setl.difference(setle)
            labelidx = [label.index(i) for i in setl]
            label = [label[i] for i in labelidx]
            sequence = [sequence[i] for i in labelidx]

            if DEBUG:
                print(f"len(label) = {len(label)}")
                print(f"len(sequence) = {len(sequence)}")
            labels.append(label)
            sequences.append(sequence)
            varsitess.append(varsites)
            labels = [[x.lower() for x in label] for label in labels]
    return [labels,sequences,varsitess,real_loci]
    
    
#====================================================================================
def read_nexus(nexus_file, ttype):
    print('\n++++++++++++ nexus file information ++++++++++++')
    with open(nexus_file,'r') as f:
        mynexus = f.readlines()
        matrix_idx=[]
        end_idx=[]
        beginsets_idx=[]
        
        label_nex = []
        sequences_nex = []
        charpartition = []
        taxpartition = []
        taxpartition_names = []

        #To extract information from nexus file based on where "matrix", "end", and "begin sets" aappear:
        for (i, line) in enumerate(mynexus):
            if 'matrix' in line:
                matrix_idx.append(i)
            if 'end' in line:
                end_idx.append(i)
            if 'begin sets;' in line:
                beginsets_idx.append(i)
                
        if DEBUG:
            print(f"\nmatrix_idx = {matrix_idx}")
            print(f"\nend_idx = {end_idx}")
            print(f"\nbeginsets_idx = {beginsets_idx}")
        
        #extract labels and sequences:
        for i in range(matrix_idx[0]+1,end_idx[0]-1 ):
            myline = mynexus[i]
            if myline=='':
                continue
            if type=='STANDARD':
                l = myline[:10]    #this assumes standard phylip format
                s = myline[11:]
            else:
                index = myline.rfind('  ')
                l = myline[:index]
                s = myline[index+1:]
            label_nex.append(l.strip().replace(' ','_'))
            sequences_nex.append(s.strip())
        if DEBUG:
            print ("\nfirst label_nex:",label_nex[0])
            print ("\nlast  label_nex:",label_nex[-1])
            print ("\nfirst sequences_nex:",sequences_nex[0])
            print ("\nlast  sequences_nex:",sequences_nex[-1])
        
        #extract charpartition loci:
        start=0
        for i in range(beginsets_idx[0]+1,end_idx[1]):
            myline = mynexus[i]
            locus_ending = myline[myline.find('-')+1:-2]
            charpartition.append([start,int(locus_ending)])
            start = int(locus_ending)
        print ("\ncharpartition = ",charpartition)
        
        #extract taxpartition loci:
        start=0
        startline= mynexus[beginsets_idx[1]+1]
        if startline.find('-') == -1:
            if DEBUG:
                print("did not find -")
            startline_idx = beginsets_idx[1]+2
        else:
            startline_idx = beginsets_idx[1]+1
        for i in range(startline_idx,end_idx[2]):
            myline = mynexus[i]
            tax_ending = myline[myline.find('-')+1:-2]
            tax_name = myline[:myline.find(':')]
            tax_name= tax_name.strip().replace(' ','')
            if DEBUG:
                print(f"tax_ending={tax_ending}")
                print(f"tax_name={tax_name}")
            taxpartition.append([start,int(tax_ending)])
            taxpartition_names.append(tax_name)
            start = int(tax_ending)
        print ("\ntaxpartition = ",taxpartition)
        print ("\ntaxpartition_names = ",taxpartition_names)
            
    #----------------labels,sequences,varsitess-------------
    labels = [label_nex]*len(charpartition)
    sequences = []
    varsitess = []
    for i,item in enumerate(charpartition):
        seqs_locus_i = []
        starting = item[0]
        ending = item[1]
        for seq in sequences_nex:
            seqs_locus_i.append(seq[starting:ending])
        sequences.append(seqs_locus_i)
        
        varsites_locus_i = [list(si) for si in seqs_locus_i if len(si)>0]
        varsites_locus_i = [len([i for i in list(set(si)) if i!='-']) for si in zip(*varsites_locus_i)]
        varsites_locus_i = [sum([vi>1 for vi in varsites_locus_i]),len(varsites_locus_i)]
        varsitess.append(varsites_locus_i)
                
    return [labels,sequences,varsitess,taxpartition,taxpartition_names]


#====================================================================================
# Use this to read labels from files for exclude or include operations and return as sets
def read_inexfiles(file_name):
    if file_name == None:
        return []
    # contains a name one by line
    with open(file_name, 'r') as f:
        return set([line.rstrip() for line in f])

#====================================================================================
# create kmers
def kmer_docs(label, sequence, varsites, kmerrange, options):
    gaps_type = options['gaps_type']
    kmers_type = options['kmers_type']

    docs=[]
    miss = 0
    
    for k in kmerrange:
        tokenize_k = tokenize(sequence, k, gaps_type, kmers_type)
        if len(docs)>0:
            docs =[docs[i]+tokenize_k[i] for i in range(len(tokenize_k))]
        else:
            docs = tokenize_k
            
    if DEBUG:
        taxa_names = label
        count = 0
        for doc in docs:
            count += len(doc)
        print(f"Number of all words in docs = {count}")
        print(f"len(taxa_names) = {len(taxa_names)} ")
        print(f"taxa_names = {taxa_names} ")
    return docs



#====================================================================================
def merge_documents(options, label, docs,taxpartition,taxpartition_names):
    datainput = options['datainput']
    merging = options['merging']
    if datainput=="folder":
        if merging:
            taxa_names = list(dict.fromkeys([x[:merging] for x in label ]))  #remove duplicates from a list,
            #while preserving order using"dict.fromkeys" insted of "list(set)"
            if DEBUG:
                print(f"distinct taxa_names = {taxa_names} ")
                
            taxa_names_indices=[]
            for item in taxa_names:
                indices= [label.index(i) for i in label if item in i]
                taxa_names_indices.append(indices)
                taxa_names_indices_len=[len(i) for i in taxa_names_indices]
                if DEBUG:
                    print(f"taxa_names_indices = {taxa_names_indices}")
                    print(f"taxa_names_indices_len = {taxa_names_indices_len}")
                docs_merged=[]
                for indices in taxa_names_indices:
                    doc_merged = list(itertools.chain.from_iterable([docs[i] for i in indices] ))
                    docs_merged.append(doc_merged)

                #DEBUG: count number of all words in docs
                if DEBUG:
                    print(f"len(docs_merged) = {len(docs_merged)} ")
                    count = 0
                    for doc in docs_merged:
                        count += len(doc)
                    print(f"Number of all words in docs_merged = {count}")
            return taxa_names, docs_merged
        else:
            taxa_names = label
            return taxa_names, docs
            
    elif datainput=="nexus_file":
        taxa_names = taxpartition_names
        if DEBUG:
            print(f"taxa_names = {taxa_names} ")
        
        docs_merged=[]
        for i,item in enumerate(taxpartition):
            doc_merged = list(itertools.chain.from_iterable(docs[item[0]:item[1]]))
            docs_merged.append(doc_merged)
        if DEBUG:
            print(f"len(docs_merged) = {len(docs_merged)} ")
            #count number of all words in docs
            count = 0
            for doc in docs_merged:
                count += len(doc)
            if DEBUG:
                print(f"Number of all words in docs_merged = {count}")
        docs = docs_merged
        return taxa_names, docs
        
        
#====================================================================================
def training(taxa_names, docs, options):
    miss = 0
    dictionary = Dictionary(docs)
    chunksize = options['chunksize']
    iterations = options['iterations']
    num_topics = options['num_topics']
    passes = options['passes']
    eval_every = options['eval_every']
    update_every = options['update_every']
    alpha = options['alpha']
    eta = options['eta']
    
    #--------------- Filtering: remove rare and common tokens ---------------
    dictionary.filter_extremes(no_below=2, no_above=0.5)
    if DEBUG:
        print(f"\nDictionary after filtering:\n{dictionary}")
        
    if len(dictionary)<1:
        if DEBUG:
            print("ZERO SIZED DICTIONARY")
        miss += 1
        return [None,None,miss]

    #---------------------- Vectorize data: Bag-of-words --------------------
    corpus = [dictionary.doc2bow(doc) for doc in docs]
    if DEBUG:
        print(f'\nNumber of documents: {len(corpus)}')
        print(f'Number of unique tokens: {len(dictionary)}')
        
    #------------------------- LDA Model: Training --------------------------
    # Make a index to word dictionary.
    temp = dictionary[0]
    id2word = dictionary.id2token
    
    model = LdaModel(corpus=corpus, id2word=id2word, chunksize=chunksize, \
                     alpha=alpha, eta=eta, \
                     iterations=iterations, num_topics=num_topics, \
                     passes=passes, eval_every=eval_every, minimum_probability=0, update_every=update_every )
    #--------------------- Print topics/words frequencies -------------------
    if DEBUG:
        for idx, topic in model.print_topics(-1):
            print("Topic: {} \nWords: {}\n".format(idx, topic))
                
    #---------------- Assigning the topics to the documents -----------------
    docs_tuples = model[corpus]
                
    #--------------------- topics list for current locus --------------------
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
    return [topics,taxa_names,miss]


#====================================================================================
def read_one_data(current, folder, locus,options):
    if DEBUG:
        print(f"current = {current}")
    prefix = options['prefix']
    suffix = options['suffix']
    ttype = options['ttype']
    filetype = options['filetype']
    include_names = options['include_names']
    exclude_names = options['exclude_names']

    labels=[]
    sequences=[]
    varsitess=[]
    seti = set(include_names)
    sete = set(exclude_names)
    real_loci = 0
    i = locus
    #Extract labels and sequences:
    locus_i = prefix+str(i)+suffix
    locus = os.path.join(current,folder,locus_i)        
    label,sequence,varsites = phylip.readData(locus, ttype)
    setl = set(label)
    if seti==set():
        seti = setl
    iok = (setl.intersection(seti) == seti)
    if iok:
        real_loci += 1
        setle = setl.intersection(sete)
        setl  = setl.difference(setle)
        labelidx = [label.index(i) for i in setl]
        label = [label[i] for i in labelidx]
        sequence = [sequence[i] for i in labelidx]
        if DEBUG:
            print(f"len(label) = {len(label)}")
            print(f"len(sequence) = {len(sequence)}")
        label = [x.lower() for x in label]
        return [label,sequence,varsites,real_loci] # a single locus
    else:
        return [None,None,None,0]
    
#====================================================================================
def use_options(current, folder, gaps_type, kmers_type, bootstrap, nbootstrap, datainput, merging, chunksize, iterations,num_topics,passes, eval_every, update_every, alpha, eta, prefix, suffix, ttype, filetype, include_names, exclude_names):
    options={}
    options['current'] = current
    options['folder'] = folder
    options['gaps_type']=gaps_type
    options['kmers_type'] = kmers_type
    options['bootstrap'] = bootstrap
    options['nbootstrap'] = nbootstrap
    options['datainput'] = datainput
    options['merging'] = merging
    options['chunksize'] = chunksize
    options['iterations']=iterations
    options['num_topics'] = num_topics
    options['passes'] = passes
    options['eval_every'] = eval_every
    options['update_every'] = update_every
    options['alpha'] = alpha
    options['eta'] = eta
    options['prefix'] = prefix
    options['suffix'] = suffix
    options['ttype'] = ttype
    options['filetype'] = filetype
    options['include_names'] = include_names
    options['exclude_names'] = exclude_names
    return options


#====================================================================================
def process_locus(i, label, sequence, varsites,taxpartition,taxpartition_names, kmerrange, options):
        
    datainput = options['datainput']
    nbootstrap = options['nbootstrap']
    bootstrap = options['bootstrap']
    gaps_type = options['gaps_type']

    #print(i,end=' ', file=sys.stderr, flush=True)
    if datainput=="folder":
        if label == []: #assume we need to read the data for each locus now and have not done that beforehand
            label, sequence, varsites, real_locus  = read_one_data(options['current'], options['folder'], i, options)
            
            
    if bootstrap == 'seq' and nbootstrap>0:
        sequence = bootstrap_sequences_one_locus(sequence)
        
 
    if gaps_type == 'rm_row':
        sequence = [seq.replace('-', '') for seq in sequence]
    elif gaps_type == 'rm_col':
        # String List to Column Character Matrix Using zip() + map()
        seq_matrix = np.transpose(list(map(list, zip(*sequence))))
        #remove columns if the column has an element '-'
        seq_matrix_cleaned = np.array([col for col in seq_matrix.T if '-' not in col]).T
        #Convert List of lists to list of Strings again
        sequence = [''.join(ele) for ele in seq_matrix_cleaned]
        
    if label != None:    
        docs = kmer_docs(label, sequence, varsites, kmerrange, options)
    
        if bootstrap == 'kmer' and nbootstrap>0:
            docs = [np.random.choice(docs[i],size=len(docs[i]),replace=True) for i in range(len(docs)) ]
    

        taxa_names, docs = merge_documents(options, label, docs, taxpartition,taxpartition_names)
        [topics,taxa_names,miss] = training(taxa_names, docs, options)
        return [topics,taxa_names,miss,real_locus]
    else:
        return [None, None, 0, 0]


#====================================================================================
def topicmodeling(options):
    global num_loci
    datainput = options['datainput']
    
    if datainput=="folder":
        kmerrange = range(kmer_range_list[0],kmer_range_list[1],kmer_range_list[2])
        args = [(i, [],[], [], [], [], kmerrange, options) for i in range(num_loci)]
    
    elif datainput=="nexus_file":
        loci = read_nexus(nexus_file, ttype)
        labels,sequences,varsitess,taxpartition,taxpartition_names = loci
        kmerrange = range(kmer_range_list[0],kmer_range_list[1],kmer_range_list[2])
        args = [(i, labels[i], sequences[i], varsitess[i],taxpartition,taxpartition_names, kmerrange, options) for i in range(num_loci)]


    mypool = multiprocessing.cpu_count() -1
    print("Number of cores to use:",mypool)
    with Pool(mypool) as p, tqdm(total=num_loci) as pbar:
        res = [p.apply_async(
            process_locus, args=args[i], callback=lambda _: pbar.update(1)) for i in range(num_loci)]
        results = [r.get() for r in res]

    
    topics_loci = []
    taxa_names_loci = []
        
    topics_loci,taxa_names_loci,miss, real_loci = zip(*results)
    topics_loci = [t for t in topics_loci if t!=None]
    taxa_names_loci  = [t for t in taxa_names_loci if t!=None]
    num_loci = np.sum(real_loci)
    miss = np.sum(miss)
    print("Missed loci out of:",miss,num_loci)
    if DEBUG:
        print(f"taxa_names_loci = {taxa_names_loci}")
    taxa_names_loci_sets = [set(l) for l in  taxa_names_loci]
    common_taxa_names = list(set.intersection(*taxa_names_loci_sets))
    uncommon_taxa_names  =  list(set.union(*taxa_names_loci_sets) - set.intersection(*taxa_names_loci_sets))
    alltaxa = list(set.union(*taxa_names_loci_sets))
    if DEBUG:
        print("ALL:", list(alltaxa))
        print("COMMON:",common_taxa_names)
        print("NOT COMMON:", uncommon_taxa_names)
            
    common_taxa_names_loci_indx = [[l.index(c) for c in common_taxa_names] for l in taxa_names_loci]
    uncommon_taxa_names_loci_indx = [[uncommoncheck(l,c) for c in uncommon_taxa_names] for l in taxa_names_loci]
    if DEBUG:
        print(f"common_taxa_names_loci_indx = {common_taxa_names_loci_indx}")
        
    common_topics_loci = [[l[i] for i in common_taxa_names_loci_indx[num]] for num,l in enumerate(topics_loci)]
    uncommon_topics_loci = [[uncommonresolve(l,i, num_topics) for i in uncommon_taxa_names_loci_indx[num]] for num,l in enumerate(topics_loci)]
        
    commonplus_taxa_names = common_taxa_names + uncommon_taxa_names
    commonplus_topics_loci = [l+u for l,u in zip(common_topics_loci, uncommon_topics_loci) ]
    if DEBUG:
        print(f"commonplus_taxa_names = {commonplus_taxa_names}")
        print(f"commonplus_topics_loci = {commonplus_topics_loci}")

    if force:
        return commonplus_taxa_names, commonplus_topics_loci,miss
    else:
        return common_taxa_names, common_topics_loci,miss


#====================================================================================
def uncommoncheck(l,c):
    if c in l:
        return l.index(c)
    else:
        return -1

def uncommonresolve(l, i, numtopics):
    if i == -1:
        return [1.0/numtopics for i in range(numtopics)]
    else:
        return l[i]
#====================================================================================
def infile(topics_loci, taxa_names, num_loci):
    if DEBUG:
        print(f"num_loci = {num_loci}")
    taxa_names_limited = [x[:10] for x in taxa_names]      #CONTML accepts population names lass than 10 letters
    num_pop= len(topics_loci)
    with  open("infile", "w") as f:
        f.write('     {}    {}\n'.format(num_pop, num_loci))
        f.write('{} '.format(num_topics)*num_loci)
        f.write('\n')
        for i in range(num_pop):
            myname = taxa_names[i]
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
    print("Contml is running...")
    os.system(f'cat contmlinput | {PROGRAMPATH}contml2 -> contml2.log')
#    os.system(f'cat contmlinput | {PROGRAMPATH}contml')
    
    #read the outtree file
    with open('outtree', 'r') as f:
        tree = f.read().replace('\n', ' ')
    return tree

#====================================================================================
def simulation(current, folder, options):
    diverge_time = float(sim_diverge)       #diverge_time=0.0/0.01/0.05/0.1/0.2
    tree_newick = '((Arb:1,Fra:1):1,Kre:1,Rom:1);'    #We need just topology to compare
    tns = dendropy.TaxonNamespace()
    true_tree = dendropy.Tree.get(data=tree_newick,schema="newick",taxon_namespace=tns)
    true_tree.encode_bipartitions()
        
    simfiles_folder = os.path.join(current,folder)
    if DEBUG:
        print(f"simfiles_folder = {simfiles_folder}")
    current = simfiles_folder
    siminfiles_trees=[]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    num_iter =[2,5,10,20,50,100]
    sims_num = list(range(0,100))
    count_equaltrue=[]
    for num in num_iter:
        num_loci = num
        distances =[]
        include_names = []
        exclude_names = []
        for sim_num in sims_num:
            siminfile_sl ='siminfile_'+str(sim_num)+'_100_'+str(diverge_time)+'_100'
            folder = os.path.join(simfiles_folder,siminfile_sl)
            if DEBUG:
                print(f"siminfile_sl = {siminfile_sl}")
                print(f"folder = {folder}")
            #sys.exit()
            options['current'] = current
            options['folder'] = folder
            if DEBUG:
                print(f"options['current'] = {options['current']}")
                print(f"\noptions = {options}")
            taxa_names, topics_loci, miss = topicmodeling(options)
            #for each locus remove last column from topic matrix
            topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
            
            
            #concatenation of the topics for all loci
            topics_loci_concatenated = topics_loci_missingLast[0]
            for i in range(1,len(topics_loci_missingLast)):
                topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]

            #generate infile
            infile(topics_loci_concatenated, taxa_names, num_loci)
            
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
def single_run(show=False,options={}):
    global num_loci

    taxa_names, topics_loci, miss = topicmodeling(options)
    #for each locus remove last column from topic matrix
    topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
        
    #concatenation of the topics for all loci
    topics_loci_concatenated = topics_loci_missingLast[0]
    for i in range(1,len(topics_loci_missingLast)):
        topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
        if DEBUG:
            print(f'topics_loci_concatenated =\n{topics_loci_concatenated}')
        
    #generate infile
    infile(topics_loci_concatenated, taxa_names, num_loci-miss)
    
    #run CONTML
    ourtree = run_contml(infile)
#    print(f"ourtree = {ourtree}")
    if show:
        #Figtree
        os.system(f"{PROGRAMPATH}figtree outtree")
    print('Effective loci use:', num_loci)
    return ourtree

        
#====================================================================================
def bootstrap_sequences(sequences):     #All loci
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
def bootstrap_sequence_one_locus(sequence):     #one locus
    num = len(sequence)
    nsites = len(sequence[0])
    pick = np.random.randint(0,nsites,nsites)
    newsequence=[]
    for n in range(num):
        newseq = "".join([locus[n][i] for i in  pick])
        newsequence.append(newseq)
    return newsequence
            


#====================================================================================
def bootstrap_run(bootstrap, options):
    global num_loci
    count_boot = 1
    outtrees=[]

    nbootstrap = options['nbootstrap']
    bootstrap = options['bootstrap']
#    print(f"bootstrapping based on '{bootstrap}'")
    
    for bi in range(nbootstrap):
        print(f"\nBootstrapping number '{bi}'")
        count_boot += 1
        taxa_names, topics_loci, miss = topicmodeling(options)
        
        #for each locus remove last column from topic matrix
        topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
        
        #concatenation of the topics for all loci
        topics_loci_concatenated = topics_loci_missingLast[0]
        for i in range(1,len(topics_loci_missingLast)):
            topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
            
        #generate infile
        infile(topics_loci_concatenated, taxa_names, num_loci-miss)
        
        #run CONTML
        ourtree = run_contml(infile)
        print(f"Bootstrap tree #{bi} = {ourtree}")
        outtrees.append(ourtree)
        
    return outtrees

#====================================================================================
if __name__ == "__main__":
    current = os.getcwd()

    if DEBUG:
        #To print progress of the training procedure on screen
        import logging
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.NOTSET)
    
    prefix = 'locus'
    suffix = '.txt'
    args = myparser() # parses the commandline arguments
    gaps_type = args.gaps_type
    merging = args.merging
    sim_diverge = args.sim_diverge
    kmers_type = args.kmers_type
    num_loci = args.num_loci
    folder = args.folder
    nexus_file = args.nexus_file
    bootstrap_type = args.bootstrap_type
    include_file = args.include_file
    exclude_file = args.exclude_file
    force = args.force
    showtree = args.showtree
    #gensim LDA arguments:
    num_topics = args.num_topics
    iterations = args.iterations
    passes = args.passes
    chunksize = args.chunksize
    eval_every = args.eval_every
    update_every = args.update_every
    alpha = args.alpha
    eta = args.eta
    if alpha.isdigit():
        alpha = int(alpha)
    if eta.isdigit():
        eta = int(eta)
                
    if folder:
        datainput = "folder"
    elif nexus_file:
        datainput = "nexus_file"
    print(f"\ndatainput : {datainput}")

    
    print('\n+++++++++++++++ LDA arguments +++++++++++++++++')
    print(f"num_topics = {num_topics}")
    print(f"iterations = {iterations}")
    print(f"passes = {passes}")
    print(f"chunksize = {chunksize}")
    print(f"eval_every = {eval_every}")
    print(f"update_every = {update_every}")
    
    print(f"type(alpha) = {type(alpha)}")
    print(f"type(eta) = {type(eta)}")
    
    print(f"alpha = {alpha}")
    print(f"eta = {eta}")
                                
    
    if DEBUG:
        print(f"\nnum_topics={num_topics}")

    nbootstrap = args.num_bootstrap
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

    include_names = read_inexfiles(include_file)
    exclude_names = read_inexfiles(exclude_file)
    
    options = use_options(current, folder, gaps_type, kmers_type, bootstrap, nbootstrap, datainput, merging, chunksize, iterations,num_topics,passes, eval_every, update_every, alpha, eta, prefix, suffix, ttype, filetype, include_names, exclude_names)
    
    #After cloning the repository, in topiccontml.py modify the PROGRAMPATH to the path that FigTree and CONTML are installed
    PROGRAMPATH = '~/bin/'
    
    # generates a random number seed for jumble in contml
    RANDOMSEED  = np.random.randint(10000,2**16,size=1)[0]
    if RANDOMSEED % 2 == 0:
        RANDOMSEED += 1
     
    start = time.time()
    #=================== if simulation analysis ======================
    if sim_diverge is not None:
        simulation(current, folder, options)
    #===========  General: if NOT simulation analysis  ===============
    elif bootstrap!= 'none':
        print(f"\nBootstrapping based on '{bootstrap}':\n")
        ourtree = single_run(False,options)
        print(f"TopicContml tree = {ourtree}\n")
        with open('best.tre','w') as btrees:
            btrees.write(ourtree+'\n')
            
        outtrees = bootstrap_run(bootstrap, options)
        with open('bootstrap_replicates.tre','w') as btrees:
            for tr in outtrees:
                btrees.write(tr+'\n')
                
        print(f"\nBootstrap trees using SumTrees")
        os.system(f"sumtrees.py --set-edges=support --decimals=0 --percentages --output-tree-filepath=bootstrap_target_best.tre --target=best.tre bootstrap_replicates.tre -q")
        os.system(f"sumtrees.py --set-edges=support --decimals=0 --percentages --output-tree-filepath=bootstrap_majority.tre bootstrap_replicates.tre -q")
        
        print(f"\n> Bootstrap trees are written into files 'bootstrap_target_best.tre' and 'bootstrap_majority.tre'")
        print(f"> Bootstrap replicates are written into file 'bootstrap_replicates.tre'")
        print(f"> TopicContml tree is written into file 'best.tre'")
    else:
        single_run(showtree,options)
    end = time.time()
    print(f"\n> Elapsed time = {end - start}\n")
    #variable not defined    print(f"\n> Tokenize_time = {tokenize_time}")
    citations(options)
