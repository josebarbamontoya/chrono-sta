# Chronological Supertree Algorithm (Chrono-STA)

Chrono-STA is a  method that incorporates the phylogenetic time dimension to construct supertrees, making it robust even when species overlap between phylogenies is low. For a more detailed explanation of the method, please read our manuscript (DOI: doi.org/...).

## Installation

Chrono-STA is suported on Unix, macOS, and Windows.

To download Chrono-STA, use git, and clone with `--recursive`. In a terminal or command prompt session, type:
```
git clone --recursive https://github.com/josebarbamontoya/chrono-sta
```

Use of Chrono-STA requires Python 3 and the following Python and R pachages. Install them before using the `chronosta.py` script.

To download and install Python 3, follow the instructions provided on the official Python website:
	https://www.python.org/downloads/

Install the required Python 3 packages. In a terminal or command prompt session, type:

	pip3 install biopython
	pip3 install pandas
	pip3 install numpy
	pip3 install scipy
	pip3 install matplotlib
	pip3 install rpy2
	pip3 install ete3 

Install the required R packages. In an R session, type:

	install.packages("ape")
	install.packages("phytools")

## Execution

**Unix and macOS**

1.	In a terminal, create a directory (e.g., `1000_genetree_speciestree`) and copy the sub-timetrees in Newick format with the extension `.nwk`, along with the `chronosta.py` script located inside the `code` directory:
```
mkdir /Users/barba/chrono-sta/examples/1000_genetree_speciestree
cp *.nwk /Users/barba/chrono-sta/examples/1000_genetree_speciestree
cp chronosta.py /Users/barba/chrono-sta/examples/1000_genetree_speciestree
```

2.	Type `cd` followed by the path to the created directory to change the working directory to the folder that contains the `chronosta.py` script:	
```
cd /Users/barba/examples/chrono-sta/1000_genetree_speciestree
```

3.	To execute the Chrono-STA analysis, type:
```
python3 chronosta.py
```

**Windows**

1.	In a command prompt, create a directory (e.g., `1000_genetree_speciestree`) and copy the sub-timetrees in Newick format with the extension `.nwk`, along with the `chronosta.py` script located inside the `code` directory:
```
mkdir C:\Users\barba\chrono-sta\examples\1000_genetree_speciestree
copy *.nwk C:\Users\barba\chrono-sta\examples\1000_genetree_speciestree
copy chronosta.py C:\Users\barba\chrono-sta\examples\1000_genetree_speciestree
```

2.	Type `cd` followed by the path to the created directory to change the working directory to the folder that contains the `chronosta.py` script:	
```
cd C:\Users\barba\chrono-sta\examples\1000_genetree_speciestree
```

3.	To execute the Chrono-STA analysis, type:
```
python3 chronosta.py
```

## Output files

- `_pd_matrix.csv` suffix denotes files containing the pairwise distance matrix computed for each sub-timetree.
- `combo_matrix.csv` contains the cumulative average pairwise distance matrix computed from sub-timetrees. 
- `clusters_and_pairwise_distances_list.csv` contains information about the clusters inferred during the construction of the supertree, along with their pairwise distances. 
- `final_pairwise_distance_matrix.csv` contains the final pairwise distance matrix computed from the list of clusters and pairwise distances.
- `supertree_from_final_pairwise_distance_matrix.newick` contains the supertree constructed from the final pairwise distance matrix in Newick format. This supertree is ultrametricized using non-negative least squares, ensuring uniform branch lengths from root to tips by minimizing the sum-of-squares distance between input and output timetrees
- `chronosta_supertimetree_pairwise_distance_matrix.csv` contains the pairwise distance matrix computed from the Chrono-STA supertimetree.
- **`chronosta_supertimetree.newick`** contains the Chrono-STA supertimetree in Newick format.
- `chronosta_out.log` contains the terminal/command prompt output produced by Chrono-STA.

---
We hope you find this repository useful. For comments and questions please e-mail jbarba@amnh.org.
