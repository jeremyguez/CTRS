# CTRS

Code for research paper in Genetics: "Cultural transmission of reproductive success impacts genomic diversity, coalescent tree topologies and demographic inferences"


## Files

`CalcIPrimeV2_linux64.exe` is used to compute the Ib index from Brandenburg et al., 2012 (input: newick trees, output: imbalance data). File obtained from J.-T. Brandenburg.

`format_newick.py` is a script written by J.-T. Brandenburg to convert tree sequences to newick format for the .exe input (used in Brandenburg et al., 2012).

`sumstats.py` contains functions for computing sumstats (using tskit, Kelleher et al., 2016) and for dadi inference (Gutenkunst et al., 2009). 

`simupop.slim` contains the SLiM code for CTRS and demographic modelisation.


## Pipeline summary

- The needed parameters for the SLiM simulation were written in a JSON file and loaded in a python script that automatically wrote them into the SLiM code. The python script then launched the SLiM code using the method `os.popen()`.

- The outputs from the SLiM simulations are tree sequences. Those were loaded in a python script. Recapitation and sampling were done, and mutations added to the trees (using tskit, Kelleher et al., 2016).

- Several summary statistics and dadi inferences were computed and the data written as .csv files. 

- An R script loaded the final data, and `ggplot2` was used to plot the papers' figures.
