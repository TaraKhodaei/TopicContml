<div align="center"><img src="https://github.com/TaraKhodaei/TPContml/blob/main/images/workflow.jpg" width="350"/></div>

Python package **TopicContml** uses $k$-mers and probabilistic topic modeling, an unsupervised machine learning approach based on natural language processing, to construct evolutionary relationships among species from unaligned DNA sequences.


# $\color{purple}{\textsf{Usage}}$
    pathtrees.py [-h] [-e] [-m] [-gt GAPS_TYPE] [-kt KMER_TYPE] [-kr KMER_RANGE] [-n NUM_LOCI] 
                        

# $\color{purple}{\textsf{Arguments}}$

**-h, --help**
> show this help message and exi  

<br/>

**-e, --extended**
> If the phylip dataset is in the extended format, use this. 
 
> <br/>

**-m, --merging**
> It merges sequences with the same population

<br/>

**-gt GAPS_TYPE, --gaps_type GAPS_TYPE**
> String "rm_row": removes gaps(-) in each sequence by the row. String "rm_col": romoves the column if there is at least one gap(-) in that column. Otherwise, it does not make changes in sequences.

<br/>

**-kt KMER_TYPE, --kmer_type KMER_TYPE**
> default "not_overlap": extract kmers without overlapping. String "not_overlap": extract kmers with overlapping.

<br/>

**-kr KMER_RANGE, --kmer_range KMER_RANGE**
> range of kmers extraction. default range is 2,10,2

<br/>

**-n NUM_LOCI, --num_loci NUM_LOCI**
> number of loci

<br/>

# $\color{purple}{\textsf{Dataset}}$
The dataset of sequences should be in the directory in a folder called **"loci"**, and inside that each locus dataset is included in separate folders called **"loci0", "loci1", ...**

<br/>

> **Experiment 1.1**.<br/>
> > **loci_birds:** The bird sequences are collected from 14 loci and 9 different locations.
> ```
> python TopicPContml.py -e -m -gt rm_row -n 19
> ```
> <div align="center"><img src="https://github.com/TaraKhodaei/TPContml/blob/main/images/experiment1_birds.jpg" width="300"/></div>





1 locus, no merging:
python TopicPContml.py -e -gt rm_row -n 1
