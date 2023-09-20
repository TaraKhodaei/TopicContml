<div align="center"><img src="images/workflow_new.jpg" width="600"/></div>

Python package **TopicContml** uses $k$-mers and probabilistic topic modeling, an unsupervised machine learning approach based on natural language processing, to construct evolutionary relationships among multilocus species from unaligned DNA sequences.$$\color{purple}{\textsf{Usage}}$$


# $$\color{purple}{\textsf{Usage}}$$
    topiccontml.py [-h] [-e] [-m MERGING] [-gt GAPS_TYPE] [-kt KMER_TYPE] [-kr KMER_RANGE]
                   [-nl NUM_LOCI] [-nt NUM_TOPICS] [-f FOLDER] [-sd SIM_DIVERGE] [-nb BOOTSTRAP] [-bt BOOTSTRAP_TYPE]
                        

# $\color{purple}{\textsf{Arguments}}$

**-h, --help**
> show this help message and exi  

<br/>

**-e, --extended**
> If the phylip dataset is in the extended format, use this. 
 
> <br/>

**-m MERGING, --merging MERGING**
> Merge sequences that start with the same number of MERGING letters.

<br/>

**-gt GAPS_TYPE, --gaps_type GAPS_TYPE**
> String "rm_row": removes gaps(-) in each sequence by the row. String "rm_col": romoves the column if there is at least one gap(-) in that column. Otherwise, it does not make changes in sequences.

<br/>

**-kt KMER_TYPE, --kmer_type KMER_TYPE**
> default "not_overlap": extract kmers without overlapping. String "overlap": extract kmers with overlapping.

<br/>

**-kr KMER_RANGE, --kmer_range KMER_RANGE**
> range of kmers extraction (default range is 2,10,2, which means it creates words with lengths 2, 4, 6, and 8).

<br/>

**-nl NUM_LOCI, --num_loci NUM_LOCI**
> number of loci

<br/>

**-nt NUM_TOPICS, --num_topics NUM_TOPICS**
> number of topics. Defult value is 5 topics.

<br/>

**-f FOLDER, --folder FOLDER**
> the folder that contains the data (loci in separate text files called "loci0.txt", "loci1.txt", ...). The defult is "loci"

<br/>

**-sd SIM_DIVERGE, --siminfile_diverge_time SIM_DIVERGE**
> To do siminfile analysis for the folder with the given float number of sim_diverge

<br/>

**-nb BOOTSTRAP, --bootstrap BOOTSTRAP**
> number of bootstrap replicates

<br/>

**-bt BOOTSTRAP_TYPE, --bootstrap_type BOOTSTRAP_TYPE**
> default "kmer": do the bootstrap by randomly choosing  x kmers in each document of x kmers. String "seq": do the bootstrap by randomly choosing  x columns  of aligned sequences with the same length of x ("seq" works only in the case the sequences have the same lengths)

<br/>

# $\color{purple}{\textsf{Requirements}}$
* The following **packages** are required: <br/>
1. `gensim`:
    ```
    pip install gensim
    ```
2. `CONTML` You will need to compile a customized version of `CONTML`, we call it `CONTML2` because if you want to run bootstrap the standard contml may fail if two individuals have the same frequencies, the custom version allows for that, it also uses a default of 15 characters for the individual names. We suggest that you create bin directory in your homedirectory and place the binaries there. The full Phylip version is here: <a html="https://evolution.genetics.washington.edu/phylip.html">https://evolution.genetics.washington.edu/phylip.html</a>

    ```
    #use this commandline snippet to compile the custom version of contml
    cd phylip-part-3.69/src
    make contml
    # you may need to use this once: mkdir -p ~/bin
    cp contml ~/bin/contml2
    ```
3. `FigTree`, download Figtree from here <a html="https://github.com/rambaut/figtree/releases">https://github.com/rambaut/figtree/releases</a>
* After cloning the repository, in `topiccontml.py` modify the `PROGRAMPATH` to the path that FigTree and CONTML are installed.
* The  **dataset** of sequences should be in the directory in a folder, all loci are in separate text files called **"loci0.txt", "loci1.txt", ...** that follow the Phylip syntax.


<br/>


# $\color{purple}{\textsf{Application to Real Data}}$

## $\color{purple}{\textsf{1. Birds Dataset}}$

> **Experiment 1.1**.<br/>
> > **loci_birds:** The bird sequences are collected from 14 loci and 9 different locations. For each locus, the length of each sequence varies from 288 to 418 base pairs, and the number of sequences varies from 78 to 92 individuals.  In each locus, we merge the words from the same location (using 3 first letters) and then apply LDA.
> ```
> python topiccontml.py -f loci_birds -m 3 -gt rm_row -nl 14
> ```
> <div align="center"><img src="images/experiment1_birds.jpg" width="300"/></div>

## $\color{purple}{\textsf{2. Simulated Dataset}}$
> **Experiment 2.1**.<br/>
> ```
> python topiccontml.py -f sim_100_0.0_100 -sd 0.0 -m 3 -gt rm_row -nl 100
> python topiccontml.py -f sim_100_0.01_100 -sd 0.01 -m 3 -gt rm_row -nl 100
> python topiccontml.py -f sim_100_0.05_100 -sd 0.05 -m 3 -gt rm_row -nl 100
> python topiccontml.py -f sim_100_0.1_100 -sd 0.1 -m 3 -gt rm_row -nl 100
> python topiccontml.py -f sim_100_0.2_100 -sd 0.2 -m 3 -gt rm_row -nl 100
> ```


<br/>


# $\color{purple}{\textsf{Bootstrap Analysis}}$
1. generate concatenated nexus file of common sequences in all loci called "`myfile.nex`", and a copy of dataset folder "loci_birds", called  "`loci_birds_copy`", such that loci are modified and contain just the common sequences in all loci:
    ```
    python nexus.py -e -nl 14 -s birdspecies -w myfile -f loci_birds
    ```
 2. Use "myfile.nex" as an input in `PAUP` to get the the SVDquartets bootstrap tree called "`svdq_tree`":
 > TopicContml> paup  <br/>
 > paup> ```exe myfile.nex```  <br/>
 > paup> ```svdq partition=birdspecies showScores=no seed=1234568 bootstrap nreps=100```  <br/>
 > paup> ```savetrees file=svdq_tree format=altnex```  <br/>

3. Use "loci_birds_copy" folder to get the TopicContml bootstrap tree: <br/>
> **No gaps** 
> ```
> python topiccontml.py -f loci_birds_copy  -m -gt rm_row -nl 14 -nb 100
> ```
> **With gaps** 
> ```
> python topiccontml.py -f loci_birds_copy  -m  -nl 14 -nb 100
> ```

