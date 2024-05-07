# Chrono-STA (Chronologial Supertree Algorithm)

Chrono-STA is a  method that incorporates the phylogenetic time dimension to construct supertrees, making it robust even when species overlap between phylogenies is low.

## Installation

Chrono-STA is suported on Windows, Unix, and macOS.

To download Chrono-STA, use git, and clone with `--recursive`:
```
git clone --recursive https://github.com/josebarbamontoya/chrono-sta
```

Use if Chrono-STA requires Python 3 and the following Python and R pachages. Install them in advance before using the program.

To download and install Python 3, follow the instructions provided on the official Python website:
	https://www.python.org/downloads/

Install the required Python packages. In a terminal or command prompt session, type:

	pip3 install os
	pip3 install Bio
	pip3 install Bio.Phylo.TreeConstruction 
	pip3 install pandas
	pip3 install numpy
	pip3 install scipy.spatial.distance
	pip3 install scipy.cluster 
	pip3 install matplotlib.pyplot 
	pip3 install csv
	pip3 install sys
	pip3 install rpy2

Install the required R packages. In an R session, type:

	install.packages("ape")
	install.packages("phytools")

## Execution

Create a directory and copy the partial timetrees in Newick format with the suffix `.nwk`, along with the `chronosta.py` program:

In a terminal or command prompt session, type `cd` followed by the path to the directory to change the working directory to the folder that contains the `chronosta.py` program:	
```
cd /Users/barba/chronosta_example
```

To execute the Chrono-STA analysis type:
```
python3 chronosta.py
```

## Output files

- suffix `_pd_matrix.csv` files contain the pairwise distance matrix computed from each partial timetree.
- `combo_matrix.csv` contains the cumulative average pairwise distance matrix computed from partial timetrees. 
- `clusters_and_pairwise_distances_list.csv` contains information about the clusters formed during the construction of the supertree, along with their pairwise distances. 
- `final_pairwise_distance_matrix.csv` contains the final pairwise distance matrix computed from the list of clusters and pairwise distances.
- `supertree_from_final_pairwise_distance_matrix.newick` contains the supertree constructed from the final pairwise distance matrix in Newick format.
- `chronosta_supertimetree_pairwise_distance_matrix.csv `contains the pairwise distance matrix computed from the Chrono-STA supertimetree.
- `chronosta_supertimetree.newick `contains the Chrono-STA supertimetree in Newick format. Note that supertree_from_final_pairwise_distance_matrix.newick underwent ultrametricization due to branch length rounding. This was done to ensure uniform branch lengths from the root to the tips in the timetree.


---
We hope you find this repository useful. For comments and questions please e-mail jbarba@amnh.org.
