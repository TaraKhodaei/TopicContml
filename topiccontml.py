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
import euclid
import dendropy
from dendropy.calculate import treecompare
import os
import sys
import itertools
import time
from itertools import chain
from collections import Counter

import pyLDAvis
import pyLDAvis.gensim_models as gensimvis
from gensim.models import CoherenceModel
import shutil
import math
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
    print("|     Journal of machine Learning research, 3:993--1022                           |")
    print("| Felsenstein, J. (2005). PHYLIP (Phylogeny Inference Package) version 3.6.       |")
    print("|     Distributed by the author. Department of Genome Sciences, University        |")
    print("|     of Washington, Seattle. (https://phylipweb.github.io/phylip/)               |")
    print("| Rehurek, R., and Sojka, P. (2010). Software framework for topic modelling with  |")
    print("|     large corpora. In proceedings of LREC 2010 Workshop on New Challenges       |")
    print("|     for NLP Frameworks, Valletta, Malta, pp.45--50.                             |")
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
                        help='Merge sequences that start with the same m letters [e.g. population or species labels].')
    parser.add_argument('-k','--kmers', dest='kmers',
                        nargs='?', const='2,10,2', default=None,
                        help='range of kmers extraction: [lowerbound,max+1,step]. Default is 2,10,2 which leads to k-mers of lengths: 2,4,6,8')
    parser.add_argument('-kt','--kmer_type', dest='kmers_type',
                        default='not_overlap', action='store',type=str,
                        help='default "not_overlap": extract kmers without overlapping. String "overlap": extract kmers with overlapping.')
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
    parser.add_argument('-neighbor','--neighbor', dest='useneighbor',
                        default=False, action='store_true',
                        help='uses Neighbor-Joining (neighbor) instead of contml')
    parser.add_argument('-tmap', dest='tmap',
                        default=False, action='store_true',
                        help='uses pyLDAvis to map the topics')
    parser.add_argument('-threads','--threads', dest='max_threads',
                        default=1000, action='store', type=int,
                        help='Number of cpu cores to use for locus-parallel runs, default is system max.')

    #gensim LDA arguments:
    parser.add_argument('-nt','--num_topics', dest='num_topics',
                        default=None, action='store', type=int,
                        help='Number of requested latent topics to be extracted from the training corpus. Defult value is 5 topics.')
    parser.add_argument('-cr','--coherence_range', dest='coherence_range',
                        nargs='?', const='2,20,4', default=None,
                        help='Compute coherence for various number of topics in range of [start, limit, step]. Defult is 2,20,4, which means to compute coherence for topics number of 2,6,10,14,18')
    parser.add_argument('-i','--iterations', dest='iterations',
                        default=500, action='store', type=int,
                        help='Maximum number of iterations through the corpus when inferring the topic distribution of a corpus. Defult value is 100 iterations.')
    parser.add_argument('-p','--passes', dest='passes',
                        default=50, action='store', type=int,
                        help='Number of passes through the corpus during training. Defult value is 5.')
    parser.add_argument('-cs','--chunksize', dest='chunksize',
                        default=2000, action='store', type=int,
                        help='Number of documents to be used in each training chunk. Defult value is 2000.')
    parser.add_argument('-ee','--eval_every', dest='eval_every',
                        default=1, action='store', type=int,
                        help='Log perplexity is estimated every that many updates. Defult value is 1.')
    parser.add_argument('-ue','--update_every', dest='update_every',
                        default=1, action='store', type=int,
                        help='Number of documents to be iterated through for each update. Defult value is 1.')
    parser.add_argument('-al','--alpha', dest='alpha',
                        default='auto',
                        help='a priori belief on document-topic distribution. It can be: (1) scalar for a symmetric prior over document-topic distribution, (2) 1D array of length equal to num_topics to denote an asymmetric user defined prior for each topic. (3) Alternatively default prior strings:"symmetric": a fixed symmetric prior of 1.0 / num_topics,"asymmetric": a fixed normalized asymmetric prior of 1.0 / (topic_index + sqrt(num_topics)),"auto":Learns an asymmetric prior from the corpus')
    parser.add_argument('-et','--eta', dest='eta',
                        default='auto',
                        help='a priori belief on topic-word distribution. It can be: (1) scalar for a symmetric prior over  topic-word distribution, (2) 1D array of length equal to num_words to denote an asymmetric user defined prior for each word, (3) matrix of shape (num_topics, num_words) to assign a probability for each word-topic combination. (4) Alternatively default prior strings:"symmetric": a fixed symmetric prior of 1.0 / num_topics,"auto": Learns an asymmetric prior from the corpus.')
                        
    parser.add_argument('-fb','--filter_below', dest='filter_below',
                        default=2, action='store', type=int,
                        help='Filter out tokens that appear in less than filter_below documents (absolute number) . Defult value is 2')                      
    parser.add_argument('-fa','--filter_above', dest='filter_above',
                        default=.5, action='store', type=float,
                        help='Filter out tokens that appear in more than filter_above documents (fraction of total corpus size, not absolute number). Defult value is 0.5')
    
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    return args
    

#====================================================================================
def use_options(current, folder, gaps_type, kmers, kmers_type, bootstrap, nbootstrap, datainput, merging, chunksize, iterations, num_topics, coherence_range,passes, eval_every, update_every, alpha, eta, prefix, suffix, ttype, filetype, include_names, exclude_names, tmap, filter_below, filter_above):
    options={}
    options['current'] = current
    options['folder'] = folder
    options['gaps_type']=gaps_type
    options['kmers'] = kmers
    options['kmers_type'] = kmers_type
    options['bootstrap'] = bootstrap
    options['nbootstrap'] = nbootstrap
    options['datainput'] = datainput
    options['merging'] = merging
    options['chunksize'] = chunksize
    options['iterations']=iterations
    options['num_topics'] = num_topics
    options['coherence_range'] = coherence_range
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
    options['tmap'] = tmap
    options['filter_below'] = filter_below
    options['filter_above'] = filter_above
    return options

#====================================================================================
def parse_coherence_range(coherence_range):
    try:
        start, limit, step = map(int, coherence_range.split(','))
        return list(range(start, limit, step))
    except ValueError:
        raise ValueError("coherence_range must be in the format 'start,limit,step' with integer values.")
        
#=================================================================================
def compute_coherence_values(dictionary, corpus, docs, id2word, coherencerange, options):
    chunksize = options['chunksize']
    alpha = options['alpha']
    eta = options['eta']
    iterations = options['iterations']
    passes = options['passes']
    eval_every = options['eval_every']
    update_every = options['update_every']

    """
    Compute c_v coherence for various number of topics

    Parameters:
    ----------
    dictionary : Gensim dictionary
    corpus : Gensim corpus
    docs : List of input documents
    limit : Max num of topics
    start : starting num of topics
    step : stepping through the range of topics by #step

    Returns:
    -------
    model_list : List of LDA topic models
    coherence_values : Coherence values corresponding to the LDA model with respective number of topics
    """
    coherence_values = []
    model_list = []
    for numtopics in coherencerange:
        model = LdaModel(corpus=corpus, id2word=id2word, chunksize=chunksize, \
                     alpha=alpha, eta=eta, \
                     iterations=iterations, num_topics=numtopics, \
                     passes=passes, eval_every=eval_every, minimum_probability=0, update_every=update_every )
        model_list.append(model)
        coherencemodel = CoherenceModel(model=model, texts=docs, dictionary=dictionary, coherence='u_mass')
        coherence_values.append(coherencemodel.get_coherence())
    return model_list, coherence_values, coherencerange
    
#====================================================================================
def plot_coherence_values(coherencevalues_loci, coherencerange_loci, subplots_per_row=1):
    """
    Plot coherence values for various numbers of topics with highlighted max values in each locus.

    Parameters:
    ----------
    coherencevalues_loci : List of lists of coherence values in each locus
    coherencerange_loci : List of lists of number of topics in each locus
    subplots_per_row : Number of subplots in each row
    """
    # Number of plots
    num_plots = len(coherencerange_loci)

    # Create subplots
    num_rows = (num_plots // subplots_per_row) + (num_plots % subplots_per_row > 0)
    fig, axs = plt.subplots(nrows=num_rows, ncols=subplots_per_row, figsize=(6 * subplots_per_row, 4 * num_rows))

    # Flatten the axes array for easy iteration
    axs = axs.flatten()

    # Plot each coherence value set with the corresponding range
    for i in range(num_plots):
        ax = axs[i]
        coherence_values = coherencevalues_loci[i]
        coherence_range = coherencerange_loci[i]

        ax.plot(coherence_range, coherence_values, marker='o', linestyle='-', color='orange')
        ax.set_title(f'Locus {i}', fontsize=13, fontweight='bold')
        ax.set_xlabel('Number of Topics')
        ax.set_ylabel('Coherence Value')

        # Highlight the maximum coherence value with a vertical purple line
        max_value = max(coherence_values)
        max_index = coherence_values.index(max_value)
        max_topic = coherence_range[max_index]
        ax.axvline(x=max_topic, color='green', linestyle='--')
        
        # Mark the max value point on the plot
        ax.scatter([max_topic], [max_value], color='green', s=90)
        
        # Annotate the maximum value and the corresponding topic number in purple with larger font size
        ax.annotate(f'Max: ({max_topic}, {max_value:.2f})', xy=(max_topic, max_value),
                    xytext=(max_topic, max_value - 0.15 * (max(coherence_values) - min(coherence_values))),
                    textcoords='data', color='green', ha='center', fontsize=12, fontweight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', pad=1))

    # Remove empty subplots if any
    for j in range(num_plots, len(axs)):
        fig.delaxes(axs[j])

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save and show the plot
    plt.savefig('coherencetopics_loci.pdf')
#    plt.show()
    
    
#====================================================================================
def tokenize(seq_list, k, kmers_type):
    """
    Tokenizes a list of sequences into k-mers based on the specified k-mers type.

    Parameters:
    seq_list (list of str): List of DNA sequences to be tokenized.
    k (int): Length of the k-mers.
    kmers_type (str): Type of k-mers to generate; either 'overlap' or 'not_overlap'.

    Returns:
    list of list of str: A list of documents, where each document is a list of k-mers from a sequence.
    """
    docs = []  # List to hold the tokenized documents

    # Generate overlapping k-mers
    if kmers_type == 'overlap':
        for seq in seq_list:
            doc = []
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                doc.append(kmer)
            docs.append(doc)

    # Generate non-overlapping k-mers
    elif kmers_type == 'not_overlap':
        for seq in seq_list:
            doc = [seq[i:i + k] for i in range(0, len(seq), k)]
            docs.append(doc)

    return docs
    
#====================================================================================
def kmer_docs(sequence, kmerslist, options):
    """
    Generates concatenated k-mer documents from sequences across a range of k-mer lengths.

    Parameters:
    sequence (list of str): List of DNA sequences to be tokenized.
    kmerslist (list of int): list of k-mer lengths to consider.
    options (dict): Contains options like 'kmers_type' (either 'overlap' or 'not_overlap').

    Returns:
    list of list of str: List of documents, where each document is a concatenated list of k-mers.
    """
    kmers_type = options['kmers_type']  # Get the k-mers type from options
    docs = []
    for k in kmerslist:
        tokenize_k = tokenize(sequence, k, kmers_type)  # Tokenize sequence with current k
        if len(docs) > 0:
            # Concatenate new k-mers to existing documents
            docs = [docs[i] + tokenize_k[i] for i in range(len(tokenize_k))]
        else:
            docs = tokenize_k  # Initialize docs with the first set of k-mers

    return docs  # Return the list of concatenated k-mer documents
#====================================================================================
def calculate_kprime(n, q=0.01, sigma={'G', 'C', 'A', 'T'}):
    """
    Compute the optimal k'  to minimize the probability of observing a random k-mer.
    Parameters:
    - n (int): sequence size.
    - q (float): Desired probability (default is 0.01).
    - sigma (set): Alphabet of k-mer characters (default is {'G', 'C', 'A', 'T'}).
    
    Returns:
    - int: Optimal k'.
    """
    # Calculate the size of the alphabet |Σ|
    sigma_size = len(sigma)
    
    # Apply the formula
    kprime = math.ceil(math.log(n * (1 - q) / q, sigma_size))
    
    return kprime

#====================================================================================
def kmerprime_docs(sequence, options):
    """
    Tokenize sequences into k-mers using a k' value calculated from the average sequence length.
    
    Parameters:
    - seq_list (list of str): List of sequences to be tokenized.
    - kmers_type (str): Type of k-mer generation ('overlap' for overlapping k-mers,
      'not_overlap' for non-overlapping k-mers).

    Returns:
    - docs (list of list of str): List of tokenized documents, where each document is a list of k-mers.
    - kprime (int): The k' value used for all sequences.
    """
    kmers_type = options['kmers_type']
    docs = []  # List to hold the tokenized documents
    
    # Calculate the average length of all sequences
    avg_length = int(sum(len(seq) for seq in sequence) / len(sequence))
    
    # Calculate k' using the average length
    kprime = calculate_kprime(avg_length)

    # Generate overlapping k-mers
    if kmers_type == 'overlap':
        for seq in sequence:
            doc = []
            for i in range(len(seq) - kprime + 1):
                kmer = seq[i:i + kprime]
                doc.append(kmer)
            docs.append(doc)

    # Generate non-overlapping k-mers
    elif kmers_type == 'not_overlap':
        for seq in sequence:
            doc = [seq[i:i + kprime] for i in range(0, len(seq), kprime)]
            docs.append(doc)
            
    return docs, kprime

#====================================================================================
def read_one_data(current, folder, locus, options):
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
#Read a nexus file
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
#Merge the sequences labels with given number
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
def get_full_topic_distribution(doc_bow, model, num_topics):
    """
    Function to get the full topic distribution for a document, by setting missing topics to 0 probability
    """
    topic_dist = model.get_document_topics(doc_bow, minimum_probability=0)
    # Ensure all topics are included, set missing topics to 0 probability
    topic_dict = dict(topic_dist)
    full_topic_dist = [(topic_id, topic_dict.get(topic_id, 0.0)) for topic_id in range(num_topics)]
    return full_topic_dist  
        
#====================================================================================
def training(taxa_names, docs, options):
    """
    Train an LDA model on the provided documents, optionally optimizing the number of topics using coherence.

    Parameters:
    taxa_names (list): List of taxa names corresponding to the documents.
    docs (list of lists): Tokenized documents (each document is a list of words).
    options (dict): Dictionary of training parameters, including LDA model options and filtering thresholds.

    Returns:
    list: Contains various elements such as topics, taxa names, model details, and dictionary information.
    """
    # Define ambiguous letters to be removed
    ambiguous_letters_to_remove = ['n','N']   #???
    # Remove words with specific ambiguous letters
    docs = [
        [word for word in doc if all(letter not in word for letter in ambiguous_letters_to_remove)]
        for doc in docs]
    
    miss = 0
    dictionary = Dictionary(docs)
    chunksize = options['chunksize']
    iterations = options['iterations']
    passes = options['passes']
    eval_every = options['eval_every']
    update_every = options['update_every']
    alpha = options['alpha']
    eta = options['eta']
    tmap = options['tmap']
    num_topics = options['num_topics']
    filter_below = options['filter_below']
    filter_above = options['filter_above']
    coherence_range = options['coherence_range']
    
    dictionary_before_filtering = len(dictionary)
    
    #--------------- Filtering: remove rare and common tokens ---------------
    """
    filter_extremes: dictionary.filter_extremes(no_below=5, no_above=0.5, keep_n=100000) (default)
    Filter out tokens that appear in
    (1)less than no_below documents (absolute number) or
    (2)more than no_above documents (fraction of total corpus size, not absolute number).
    After (1) and (2), keep only the first keep_n most frequent tokens (or keep all if None), the default is 100000.
    """
    dictionary.filter_extremes(no_below=filter_below, no_above=filter_above)
    dictionary_after_filtering = len(dictionary)
    
    if len(dictionary)<1:
        if DEBUG:
            print("ZERO SIZED DICTIONARY")
        miss += 1
        #return [None,None,miss]
        return [None, None,miss,0,0,0,0,0,0,0,0,0]

    #---------------------- Vectorize data: Bag-of-words --------------------
    corpus = [dictionary.doc2bow(doc) for doc in docs]
    if DEBUG:
        print(f'\nNumber of documents: {len(corpus)}')
        print(f'Number of unique tokens: {len(dictionary)}')
        
    #------------------------- LDA Model: Training --------------------------
    # Make a index to word dictionary.
    temp = dictionary[0]
    id2word = dictionary.id2token
    
    if coherence_range is not None:
       	coherencerange = parse_coherence_range(coherence_range)    #this gives a list like [2, 6, 10, 14, 18, 22, 26] for start,limit,step=2,30,4
        
        model_list, coherence_values, coherencerange = compute_coherence_values(dictionary, corpus, docs, id2word, coherencerange, options)
        
        best_coherence_index = coherence_values.index(max(coherence_values))
        num_topics = coherencerange[best_coherence_index]
        model = model_list[best_coherence_index]
    
    else:
        model = LdaModel(corpus=corpus, id2word=id2word, chunksize=chunksize, \
                        alpha=alpha, eta=eta, \
                        iterations=iterations, num_topics=num_topics, \
                        passes=passes, eval_every=eval_every, minimum_probability=0, update_every=update_every )

    # Get the full topic distributions for all documents
    docs_tuples = [get_full_topic_distribution(doc_bow, model, model.num_topics) for doc_bow in corpus]

    # Extract topics list for each document
    topics = [[y for (x, y) in doc_topic_dist] for doc_topic_dist in docs_tuples]

    
    if coherence_range:
        return [topics, taxa_names, miss, model, corpus, dictionary, num_topics, coherence_values, coherencerange, dictionary_before_filtering, dictionary_after_filtering, docs]
    else:
        return [topics, taxa_names, miss, model, corpus, dictionary, num_topics,0,0,dictionary_before_filtering, dictionary_after_filtering, docs]
        
        
#====================================================================================
def process_locus(i, label, sequence, varsites,taxpartition,taxpartition_names, options):
        
    datainput = options['datainput']
    nbootstrap = options['nbootstrap']
    bootstrap = options['bootstrap']
    gaps_type = options['gaps_type']
    coherence_range = options['coherence_range']
    kmers = options['kmers']

    if datainput=="folder":
        if label == []: #assume we need to read the data for each locus now and have not done that beforehand
            label, sequence, varsites, real_locus  = read_one_data(options['current'], options['folder'], i, options)
            
    if sequence == [] or sequence == None:
        return [None, None, 0, 0, 0, 0, 0, 0, 0,0,0,0]   #???
                
    if bootstrap == 'seq' and nbootstrap>0:
        sequence = bootstrap_sequences_one_locus(sequence)
 
    if gaps_type == 'rm_row':
        sequence = [seq.replace('-', '') for seq in sequence]
        
    
    elif gaps_type == 'rm_col':
        if DEBUG:
            #---------------------test----------------------
            # Determine the length of the longest string to check each position
            max_length = max(len(s) for s in sequence)
            
            # Create a nested list where each sublist will store the string indices with hyphens at each position
            position_hyphen_indices = []
            
            # Iterate through each character position up to the length of the longest string
            for pos in range(max_length):
                current_position_indices = []  # List to store indices of strings with a hyphen at this position
                for idx, string in enumerate(sequence):
                    # Check if the position is valid within the string length and if the character is a hyphen
                    if pos < len(string) and string[pos] == '-':
                        current_position_indices.append(idx)
                        # Append the list of indices for this position (empty if no hyphen at this position)
                        position_hyphen_indices.append(current_position_indices)
                        
                        # Filter to report only positions where there are no hyphens (empty lists)
                        positions_without_hyphen = [i for i, hyphen_indices in enumerate(position_hyphen_indices) if not hyphen_indices]
                        
                        # Return the positions without any hyphen
                        print(f"Locus {i}:\npositions_without_hyphen = {positions_without_hyphen}")
                        #---------------------test----------------------


        # String List to Column Character Matrix Using zip() + map()
        seq_matrix = np.transpose(list(map(list, zip(*sequence))))
        #remove columns if the column has an element '-'
        seq_matrix_cleaned = np.array([col for col in seq_matrix.T if '-' not in col]).T
        #Convert List of lists to list of Strings again
        sequence = [''.join(ele) for ele in seq_matrix_cleaned]
    #@@@`
    if sequence == [] or sequence == None:
        return [None, None, 0, 0, 0, 0, 0, 0, 0,0,0,0]   #???
    #@@@
    if label != None:
        start_docs = time.time()
        if kmers:
            if len(kmers)>1:    # range of kmers given by user
                kmerslist = range(kmers[0],kmers[1],kmers[2])
            elif len(kmers)==1: # one kmer given by user
                kmerslist = kmers
            docs = kmer_docs(sequence, kmerslist, options)
        else:
            docs, kprime = kmerprime_docs(sequence, options)
        
        if bootstrap == 'kmer' and nbootstrap>0:
            docs = [np.random.choice(docs[i],size=len(docs[i]),replace=True) for i in range(len(docs)) ]
                
        end_docs = time.time()
        words_time = end_docs - start_docs
#        print(f"\n> extract words run time = {end_docs - start_docs}\n")

        taxa_names, docs = merge_documents(options, label, docs, taxpartition,taxpartition_names)
        
        if coherence_range:
            [topics, taxa_names, miss, model, corpus, dictionary, num_topics, coherence_values, coherencerange, dictionary_before_filtering, dictionary_after_filtering, docs] = training(taxa_names, docs, options)

        else:
            zzzzz = training(taxa_names, docs, options)
            [topics, taxa_names, miss, model, corpus, dictionary, numtopics, dummy1, dummy2,dictionary_before_filtering, dictionary_after_filtering, docs] = zzzzz

        
        if topics == [] or topics == None:
            topics = None
        if taxa_names == [] or taxa_names == None:
            taxa_names = None
            
        
        if coherence_range:
            if kmers:
                return [topics, taxa_names, miss, real_locus, model, corpus, dictionary, num_topics, coherence_values, coherencerange, dictionary_before_filtering, dictionary_after_filtering, docs, words_time]
            else:
                return [topics, taxa_names, miss, real_locus, model, corpus, dictionary, num_topics, coherence_values, coherencerange, dictionary_before_filtering, dictionary_after_filtering, docs,  kprime, words_time]
        else:
            if kmers:
                return [topics,taxa_names,miss,real_locus, model, corpus, dictionary, dictionary_before_filtering, dictionary_after_filtering,docs, words_time]
            else:
                return [topics,taxa_names,miss,real_locus, model, corpus, dictionary, dictionary_before_filtering, dictionary_after_filtering,docs, kprime, words_time]
        
    else:
        return [None, None, 0, 0, 0, 0, 0, 0, 0,0,0,0]   #???

#====================================================================================
def topicmodeling(options):
    global num_loci
    datainput = options['datainput']
    

    if datainput=="folder":
        args = [(i, [],[], [], [], [], options) for i in range(num_loci)]
    
    elif datainput=="nexus_file":
        loci = read_nexus(nexus_file, ttype)
        labels,sequences,varsitess,taxpartition,taxpartition_names = loci
        args = [(i, labels[i], sequences[i], varsitess[i],taxpartition,taxpartition_names, options) for i in range(num_loci)]


    #mypool = multiprocessing.cpu_count() -1
    print("\nmultilocus LDA is running ...")
    print("number of cores to use:",mypool)
    
    start_lda = time.time()
    with Pool(mypool) as p, tqdm(total=num_loci) as pbar:
        res = [p.apply_async(
            process_locus, args=args[i], callback=lambda _: pbar.update(1)) for i in range(num_loci)]
        results = [r.get() for r in res]
    end_lda = time.time()
    print(f"> LDA run time = {end_lda - start_lda}")

    if coherence_range:
        if kmers:
            topics_loci, taxa_names_loci, miss, real_loci, model_loci, corpus_loci, dictionary_loci, num_topics_loci, coherencevalues_loci, coherencerange_loci, dictionary_before_filtering_loci, dictionary_after_filtering_loci, docs_loci, words_time_loci = zip(*results)
        else:
            topics_loci, taxa_names_loci, miss, real_loci, model_loci, corpus_loci, dictionary_loci, num_topics_loci, coherencevalues_loci, coherencerange_loci, dictionary_before_filtering_loci, dictionary_after_filtering_loci, docs_loci, kprime_loci, words_time_loci = zip(*results)
            kprime_loci = [t for t in kprime_loci if t!=None]
            print(f"> optimal kmer lengths across loci = {kprime_loci}")

        num_topics_loci = [t for t in num_topics_loci if t!=None]
        coherencevalues_loci = [t for t in coherencevalues_loci if t!=None]
        coherencerange_loci = [t for t in coherencerange_loci if t!=None]
        print(f"\ncoherencevalues_loci = {coherencevalues_loci}")
        print(f"coherencerange_loci = {coherencerange_loci}")
        print(f"num_topics_loci = {num_topics_loci}")
        
        
        #plot the coherence values in each locus:
        plot_coherence_values(coherencevalues_loci, coherencerange_loci, subplots_per_row=3)
        '''
        #-----------------------map topic using pyLDAvis ------------------------
        if tmap:
            folder_name = 'pyLDAvis_plots_loci'
            if os.path.exists(folder_name):
                # Empty the folder if it exists
                shutil.rmtree(folder_name)
            os.makedirs(folder_name)
            
            for i in range(len(model_loci)):
                vis_data = gensimvis.prepare(model_loci[i], corpus_loci[i], dictionary_loci[i], mds='mmds')
                visfile = os.path.join(folder_name, f'topicsmap_locus{i}.html')
                pyLDAvis.save_html(vis_data, visfile)
        '''
        
    else:
        if kmers:
            topics_loci,taxa_names_loci,miss, real_loci, model_loci, corpus_loci, dictionary_loci, dictionary_before_filtering_loci, dictionary_after_filtering_loci, docs_loci, words_time_loci = zip(*results)
        else:
            topics_loci,taxa_names_loci,miss, real_loci, model_loci, corpus_loci, dictionary_loci, dictionary_before_filtering_loci, dictionary_after_filtering_loci, docs_loci, kprime_loci, words_time_loci = zip(*results)
            kprime_loci = [t for t in kprime_loci if t!=None]
            print(f"> optimal kmer across loci = {kprime_loci}")
      
      
    print(f"> extracting words across loci time = {sum(words_time_loci)}")
    print(f"> number of words across loci before filtering = {dictionary_before_filtering_loci}")
    print(f"> number of words across loci after filtering = {dictionary_after_filtering_loci}")
    
    
    #-----------------------map topic using pyLDAvis ------------------------
    if tmap:
        folder_name = 'pyLDAvis_plots_loci'
        if os.path.exists(folder_name):
            # Empty the folder if it exists
            shutil.rmtree(folder_name)
        os.makedirs(folder_name)
        
        for i in range(len(model_loci)):
            vis_data = gensimvis.prepare(model_loci[i], corpus_loci[i], dictionary_loci[i], mds='mmds')
            visfile = os.path.join(folder_name, f'topicsmap_locus{i}.html')
            pyLDAvis.save_html(vis_data, visfile)
     #----------------------------------------------------------------------           

    topics_loci = [t for t in topics_loci if t!=None]
        
    taxa_names_loci  = [t for t in taxa_names_loci if t!=None]
    num_loci = np.sum(real_loci)
    miss = np.sum(miss)
    print("> missed loci out of:",miss,num_loci)
    if DEBUG:
        print(f"\ntaxa_names_loci = {taxa_names_loci}")
    taxa_names_loci_sets = [set(l) for l in  taxa_names_loci]
    common_taxa_names = list(set.intersection(*taxa_names_loci_sets))    #common taxa among loci
    uncommon_taxa_names  =  list(set.union(*taxa_names_loci_sets) - set.intersection(*taxa_names_loci_sets)) #uncommon taxa among loci
    alltaxa = list(set.union(*taxa_names_loci_sets))
    if DEBUG:
        print("\nALL:", list(alltaxa))
        print("\nCOMMON:",common_taxa_names)
        print("\nNOT COMMON:", uncommon_taxa_names)
     
     
    common_taxa_names_loci_indx = [[l.index(c) for c in common_taxa_names] for l in taxa_names_loci]
    uncommon_taxa_names_loci_indx = [[uncommoncheck(l,c) for c in uncommon_taxa_names] for l in taxa_names_loci]
    if DEBUG:
        print(f"\ncommon_taxa_names_loci_indx = {common_taxa_names_loci_indx}")
        
    common_topics_loci = [[l[i] for i in common_taxa_names_loci_indx[num]] for num,l in enumerate(topics_loci)]
    uncommon_topics_loci = [[uncommonresolve(l,i, num_topics) for i in uncommon_taxa_names_loci_indx[num]] for num,l in enumerate(topics_loci)]
        
    commonplus_taxa_names = common_taxa_names + uncommon_taxa_names
    commonplus_topics_loci = [l+u for l,u in zip(common_topics_loci, uncommon_topics_loci) ]
    if DEBUG:
        print(f"\ncommonplus_taxa_names = {commonplus_taxa_names}")
        print(f"\ncommonplus_topics_loci = {commonplus_topics_loci}")

    if coherence_range:
        if force:
            return commonplus_taxa_names, commonplus_topics_loci, miss, num_topics_loci
        else:
            return common_taxa_names, common_topics_loci, miss, num_topics_loci
    
    else:
        if force:
            return commonplus_taxa_names, commonplus_topics_loci, miss, num_topics
        else:
            return common_taxa_names, common_topics_loci, miss, num_topics


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
def find_and_modify_duplicates_decimal(input_list, decimal_places):
    # Dictionary to store the indexes of each element
    index_dict = {}
    # Convert input list to strings with specified decimal places for comparison
    formatted_list = [f"{x:.{decimal_places}f}" for x in input_list]
    
    # Iterate over the list with index
    for index, element in enumerate(formatted_list):
        # If the element is already in the dictionary, append the index
        if element in index_dict:
            index_dict[element].append(index)
        else:
            # Otherwise, create a new entry with the index in a list
            index_dict[element] = [index]
    modify_formatted_list = formatted_list[:]
    position = -1
    # Modify duplicates by adding decimal values
    for key, indexes in index_dict.items():
        # Skip the first occurrence, start modifying from the second one
        num = 1
        for i in range(1, len(indexes)):
            # Convert the formatted string to a list to modify the specific decimal position
            modify_key = list(key)
            new_digit = str(int(modify_key[position]) + num)
            modify_key[position] = new_digit
            # Reconstruct the modified value as a string
            modified_key = ''.join(modify_key)
            modify_formatted_list[indexes[i]] = modified_key
            num +=1
    return modify_formatted_list

#====================================================================================
def modify_infile(nested_list, decimal_places):
    # Create a deep copy to modify
    modify_nestedlist = [sublist[:] for sublist in nested_list]
    
    # Loop over each index of the sublists
    for index in range(len(modify_nestedlist[0])):
        values_count = {}  # To track occurrences of each value
        column = [sublist[index] for sublist in modify_nestedlist]
        modified_column = find_and_modify_duplicates_decimal(column, decimal_places)
        
        # Replacing the first column
        for i in range(len(modify_nestedlist)):
            modify_nestedlist[i][index] = modified_column[i]
            
    return modify_nestedlist
    
#====================================================================================
def infile_func(topics_loci, taxa_names, num_loci, numsoftopics):
    modified_topics_loci = modify_infile(topics_loci, 15)
    if DEBUG:
        print(f"num_loci = {num_loci}")
    taxa_names_limited = [x[:10] for x in taxa_names]      #CONTML accepts population names lass than 10 letters
    num_pop= len(topics_loci)
        
    with  open("infile", "w") as f:
        f.write('     {}    {}\n'.format(num_pop, num_loci))
        if coherence_range:
            f.write(' '.join(map(str, numsoftopics)) + ' ')
        else:
            f.write('{} '.format(num_topics)*num_loci)
        f.write('\n')
        for i in range(num_pop):
            myname = taxa_names[i]
            formatted_values = " ".join(f"{x}" for x in modified_topics_loci[i])
            f.write(f'{myname:<15}{formatted_values}\n')
    f.close()

#====================================================================================
def run_contml(infile):
    os.system('rm outfile outtree')
    contmljumble = RANDOMSEED
    contmltimes  = 10
    with open('contmlinput','w') as f:
        contmlinput = f'g\nj\n{contmljumble}\n{contmltimes}\ny'
        f.write(contmlinput)
        
    print("\n\nContml is running ...")
    start_contml = time.time()
    os.system(f'cat contmlinput | {PROGRAMPATH}contml2 -> contml2.log')
    end_contml = time.time()
    print(f"> Contml run time = {end_contml - start_contml}")
    
    #read the outtree file
    with open('outtree', 'r') as f:
        tree = f.read().replace('\n', ' ')
    return tree
    
#====================================================================================
def run_neighbor(infile):
    os.system('rm outfile outtree')
    print("Euclidiean distance generation and Neighbor is running...")
    outfile = INFILE
    euclid.distances(infile,outfile)
    
    os.system(f'echo "Y" | {PROGRAMPATH}neighbor -> neighbor.log')
    
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
            
            if coherence_range:
                taxa_names, topics_loci, miss, num_topics_loci = topicmodeling(options)
            else:
                taxa_names, topics_loci, miss = topicmodeling(options)
              
            #for each locus remove last column from topic matrix
            topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
            
            #concatenation of the topics for all loci
            topics_loci_concatenated = topics_loci_missingLast[0]
            for i in range(1,len(topics_loci_missingLast)):
                topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]

            #generate infile
            infile_func(topics_loci_concatenated, taxa_names, num_loci)
            if not useneighbor:
                #run CONTML
                ourtree = run_contml(infile)
            else:
                shutil.copy(infile,infile+"savecopy")
                ourtree = run_neighbor(infile)
                
            if DEBUG:
                print(f"outtree = {ourtree}")
            
            our_tree = dendropy.Tree.get(data=ourtree,schema="newick",taxon_namespace=tns)
            our_tree.encode_bipartitions()
            distance=treecompare.unweighted_robinson_foulds_distance(true_tree, our_tree, is_bipartitions_updated=True)
            distances.append(distance)
        count_equaltrue.append(distances.count(0))
            
    np.savetxt('trueAgreement_simulation{}'.format(diverge_time), count_equaltrue,  fmt='%s')
    print(f"\nResults of agreement with true tree using RF-distanc is witten to the file 'trueAgreement_simulation{diverge_time}'")

#====================================================================================
# a single analysis that shows the tree using figtree with show=True
def single_run(show=False, options={}):
    global num_loci

    taxa_names, topics_loci, miss, numsoftopics = topicmodeling(options)
                
    #for each locus remove last column from topic matrix
    topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
        
    #concatenation of the topics for all loci
    topics_loci_concatenated = topics_loci_missingLast[0]
    for i in range(1,len(topics_loci_missingLast)):
        topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
        if DEBUG:
            print(f'topics_loci_concatenated =\n{topics_loci_concatenated}')
    
    #generate infile
    infile_func(topics_loci_concatenated, taxa_names, num_loci-miss, numsoftopics)

    if not useneighbor:
        #run CONTML
        ourtree = run_contml(infile)
    else:
        shutil.copy(infile,infile+"savecopy")
        ourtree = run_neighbor(infile)
        
    if DEBUG:
        print(f"outtree = {ourtree}")
    
    if show:
        #Figtree
        os.system(f"{PROGRAMPATH}figtree outtree")
    print(F"\n> Effective loci use = {num_loci}")

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
    outtrees=[]

    nbootstrap = options['nbootstrap']
    bootstrap = options['bootstrap']
    
    for bi in range(nbootstrap):
        print(f"\nBootstrapping number '{bi}'")
        #taxa_names, topics_loci, miss = topicmodeling(options)
        taxa_names, topics_loci, miss, numsoftopics = topicmodeling(options)
        
        #for each locus remove last column from topic matrix
        topics_loci_missingLast = [[item[i][:-1] for i in range(len(item))] for item in topics_loci]
        
        #concatenation of the topics for all loci
        topics_loci_concatenated = topics_loci_missingLast[0]
        for i in range(1,len(topics_loci_missingLast)):
            topics_loci_concatenated = [a+b for a, b in zip(topics_loci_concatenated, topics_loci_missingLast[i]) ]
            
        #generate infile
        #infile_func(topics_loci_concatenated, taxa_names, num_loci-miss)
        infile_func(topics_loci_concatenated, taxa_names, num_loci-miss, numsoftopics)

        if not useneighbor:
            #run CONTML
            ourtree = run_contml(infile)
        else:
            shutil.copy(infile,infile+"savecopy")
            ourtree = run_neighbor(infile)
        
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
    tmap = args.tmap
    max_threads = args.max_threads
    useneighbor = args.useneighbor # if true use neighbor if false use contml [default]
    mypool = multiprocessing.cpu_count() -1
    if max_threads < mypool:
        mypool = max_threads
    
    #gensim LDA arguments:
    num_topics = args.num_topics
    coherence_range = args.coherence_range
    iterations = args.iterations
    passes = args.passes
    chunksize = args.chunksize
    eval_every = args.eval_every
    update_every = args.update_every
    alpha = args.alpha
    eta = args.eta
    filter_below = args.filter_below
    filter_above = args.filter_above
    kmers = args.kmers
        
    if alpha.isdigit():
        alpha = int(alpha)
    if eta.isdigit():
        eta = int(eta)
                
    if folder:
        datainput = "folder"
    elif nexus_file:
        datainput = "nexus_file"
    print(f"\ndatainput : {datainput}")
    
        
    if coherence_range is not None:
        if coherence_range == '':  # Handle the case where '-cr' is provided without a value
            coherence_range = '2,20,4'
        else:
            coherence_range = coherence_range
    else:
        num_topics = num_topics or 5

    
    print('\n============ Arguments =============')
    print(f"num_topics = {num_topics}")
    print(f"coherence_range = {coherence_range}")
    print(f"iterations = {iterations}")
    print(f"passes = {passes}")
    print(f"chunksize = {chunksize}")
    print(f"eval_every = {eval_every}")
    print(f"update_every = {update_every}")
    print(f"alpha = {alpha}")
    print(f"eta = {eta}")
    print(f"filter_below = {filter_below}")
    print(f"filter_above = {filter_above}")
    print(f"kmers_type = {kmers_type}")
    
    if gaps_type is None:
        print(f"gaps_type = {gaps_type} (gaps are NOT REMOVED)")
    else:
        print(f"gaps_type = {gaps_type} (gaps are REMOVED)")
    
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
        


    if kmers is not None:
        kmers = kmers
        kmers = list(map(int, kmers.split(',')))
        if len(kmers)>1:    # range of kmers given by user
            print("kmer length list = ",list(range(kmers[0],kmers[1],kmers[2])))
        elif len(kmers)==1: # one kmer given by user
            print("kmer length list = ", kmers)
            
    if bootstrap!= 'none':
        print("bootstrap = ", bootstrap)
    print('====================================')
    
    include_names = read_inexfiles(include_file)
    exclude_names = read_inexfiles(exclude_file)
    
    options = use_options(current, folder, gaps_type, kmers, kmers_type, bootstrap, nbootstrap, datainput, merging, chunksize, iterations,num_topics,coherence_range, passes, eval_every, update_every, alpha, eta, prefix, suffix, ttype, filetype, include_names, exclude_names, tmap, filter_below, filter_above)


    
    # sets name for data interaction with phylip programs contml2 or neighbor
    INFILE = "infile"
    infile = INFILE
    
    #After cloning the repository, in topiccontml.py modify the PROGRAMPATH to the path that FigTree and CONTML are installed
    PROGRAMPATH = '~/bin/'
    
    # generates a random number seed for jumble in contml
    RANDOMSEED  = np.random.randint(10000,2**16,size=1)[0]
    if RANDOMSEED % 2 == 0:
        RANDOMSEED += 1
     
    start_total = time.time()
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
    end_total = time.time()
    print(f"> Elapsed time = {end_total - start_total}\n")
    citations(options)
