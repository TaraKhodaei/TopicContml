<div align="center"><img src="https://github.com/TaraKhodaei/TPContml/blob/main/images/workflow.jpg" width="350"/></div>

Python package **TopicContml** uses $k$-mers and probabilistic topic modeling, an unsupervised machine learning approach based on natural language processing, to construct evolutionary relationships among species from unaligned DNA sequences.




# $\color{purple}{\textsf{Usage}}$
    pathtrees.py [-h] [-e] [-g] [-m] [-kr kmer_range]
                        [-kt kmer_type] [-n num_loci] 
                        
                        
19 loci, merging populations in each locus:
python TPContml.py -e -m -g rm_row  -kr 2,10,2 -n 19


1ocus, no merging:
python TPContml.py -e -g rm_row  -kr 2,10,2 -n 1
