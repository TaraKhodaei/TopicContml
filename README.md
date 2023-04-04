<div align="center"><img src="https://github.com/TaraKhodaei/TPContml/blob/main/images/workflow.jpg" width="350"/></div>

Python package **TopicContml** uses $k$-mers and probabilistic topic modeling, an unsupervised machine learning approach based on natural language processing, to construct evolutionary relationships among species from unaligned DNA sequences.

Python package **TPContml** enables the construction, visualization and exploration of the continuous tree landscape interior of the convex hull of given starting trees, using insights from the Billera-Holmes-Vogtmann treespace.



# $\color{purple}{\textsf{Usage}}$
    pathtrees.py [-h] [-o OUTPUTDIR] [-v] [-p PLOTFILE] [-n NUMPATHTREES]
                        [-b NUMBESTTREES] [-r NUM_RANDOM_TREES] [-g OUTGROUP]
                        [-i NUM_ITERATIONS] [-e] [-hull] [-gtp] [-nel] [-c COMPARE_TREES]
                        [-interp INTERPOLATION] [-valid]
                        STARTTREES DATAFILE
                        
                        
                        19 loci, merging populations in each locus:
python TPContml.py -e -m -g rm_row  -kr 2,10,2 -n 19


1ocus, no merging:
python TPContml.py -e -g rm_row  -kr 2,10,2 -n 1
