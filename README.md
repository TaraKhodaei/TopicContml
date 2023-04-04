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

> **Experiment 1.1**.<br/>
> ```
> python TopicPContml.py -e -m -gt rm_row -n 19
> ```
> <div align="center"><img src="https://github.com/TaraKhodaei/TPContml/blob/main/images/experiment1_birds.jpg" width="300"/></div>


19 loci, merging populations in each locus:
python TPContml.py -e -m -g rm_row  -kr 2,10,2 -n 19


1ocus, no merging:
python TPContml.py -e -g rm_row  -kr 2,10,2 -n 1
