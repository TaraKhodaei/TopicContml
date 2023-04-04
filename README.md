<div align="center"><img src="https://raw.githubusercontent.com/TaraKhodaei/TPContml/main/images/workflow.jpg" width="35"/></div>
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
