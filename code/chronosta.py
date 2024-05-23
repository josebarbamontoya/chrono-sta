#!/usr/bin/python3

####################################################################################################
##### Chronological Supertree Algorithm (Chrono-STA) ###############################################
##### Jose Barba, Jack Craig, and Sudhir Kumar #####################################################
####################################################################################################

##### Usage
#python3 chronosta.py

import os
import warnings
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import csv
import sys
import rpy2
from ete3 import Tree
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def main():
    try:
        ##### get the directory containing the script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ##### set the working directory to the directory containing the script
        os.chdir(script_dir)

        ####################################################################################################
        ##### part01 #######################################################################################
        ##### build initial/partial time distance matrices #################################################
        ####################################################################################################

        ##### get all .nwk files in the directory
        tree_files = [file for file in os.listdir('.') if file.endswith('.nwk')]

        ##### print number of tree files
        number_of_tree_files = len(tree_files)
        print(f"Number of trees in the set: {number_of_tree_files}")

        ##### designate all tip labels and make a list
        tip_labels = set()

        ##### loop through the tree files
        for file in tree_files:
            tree = Phylo.read(file, 'newick')

            ##### get the tip labels from the tree
            for clade in tree.find_clades():
                if clade.is_terminal():
                    tip_labels.add(clade.name)

        ##### convert the set to a sorted list
        tip_labels_list = sorted(list(tip_labels))

        ##### print number of unique labels
        number_of_tip_labels = len(tip_labels_list)
        print(f"Number of unique tip labels: {number_of_tip_labels}")

        ##### print the list of unique tip labels
        print("Unique tip labels:")
        print(tip_labels_list)

        ##### create an empty list to store the matrices
        matrices = []

        ##### loop through each tree file
        for tree_file in tree_files:
            #print(f"Processing tree: {tree_file}")

            ##### load the tree from a file
            tree = Phylo.read(tree_file, "newick")

            ##### compute pairwise distances
            pairwise_distances = {}
            for leaf1 in tree.get_terminals():
                for leaf2 in tree.get_terminals():
                    if leaf1 != leaf2:
                        distance = tree.distance(leaf1, leaf2)
                        pairwise_distances[(leaf1.name, leaf2.name)] = distance

            ##### create a DataFrame from pairwise_distances
            df = pd.DataFrame.from_dict(pairwise_distances, orient='index', columns=['distance'])

            ##### extract the node names from the index
            df[['sis_a', 'sis_b']] = pd.DataFrame(df.index.tolist(), index=df.index)

            ##### pivot the DataFrame to create the square matrix
            matrix = df.pivot(index='sis_a', columns='sis_b', values='distance')
            #print("square_matrix:")
            #print(matrix)

            ##### export matrix as csv file
            matrix.to_csv(f'{os.path.splitext(tree_file)[0]}_pd_matrix.csv', index=True, header=True, na_rep='NaN')

            ##### append the matrix to the list
            matrices.append(matrix)

        ####################################################################################################
        ##### part02 #######################################################################################
        ##### combine matrices and expand to match combined matrix #########################################
        ####################################################################################################

        ##### create an empty combo DataFrame with dimensions based on the length of tip_labels_list
        combo_matrix = pd.DataFrame(index=tip_labels_list, columns=tip_labels_list)

        ##### fill the DataFrame with NaN values
        combo_matrix.fillna(value=np.nan, inplace=True)

        ##### loops for expanding matrices from matrices list to match combo_matrix
        expanded_matrices = []

        for matrix in matrices:
            ##### create a new DataFrame with dimensions based on tip_labels_list
            expanded_matrix = pd.DataFrame(index=tip_labels_list, columns=tip_labels_list)

            ##### fill the expanded matrix with values from the original matrix
            for row in tip_labels_list:
                for column in tip_labels_list:
                    ##### check if the row and column labels exist in the original matrix
                    if row in matrix.index and column in matrix.columns:
                        ##### copy the value from the original matrix
                        expanded_matrix.loc[row, column] = matrix.loc[row, column]
                    else:
                        ##### fill NaN for missing values
                        expanded_matrix.loc[row, column] = np.nan

            ##### append the expanded matrix to the list
            expanded_matrices.append(expanded_matrix)

        for index, expanded_matrix in enumerate(expanded_matrices):
            ##### generate a unique file name for each expanded matrix
            file_name = f'expanded_matrix_{index}.csv'

            ##### export the expanded matrix as a csv file
            #expanded_matrix.to_csv(file_name, index=True, header=True, na_rep='NaN')

        ####################################################################################################
        ##### part03 #######################################################################################
        ##### populate expanded matrices, compute final matrix, and clusters list ##########################
        ####################################################################################################

        ##### generate an object oriented cummulative average matrix from partial matices 
        partial_matrices = expanded_matrices

        class Matrix:
            def __init__(self, data):
                self.data = data

            def fill_missing_values(self, fill_value=np.nan):
                self.data.fillna(value=fill_value, inplace=True)

            def expand_to_match(self, target_labels):
                expanded_matrix = pd.DataFrame(index=target_labels, columns=target_labels)
                expanded_matrix.fillna(value=np.nan, inplace=True)

                for row in target_labels:
                    for column in target_labels:
                        if row in self.data.index and column in self.data.columns:
                            expanded_matrix.loc[row, column] = self.data.loc[row, column]
                        else:
                            expanded_matrix.loc[row, column] = np.nan

                self.data = expanded_matrix

            def combine_matrices(self, matrices):
                combo_matrix = pd.DataFrame(0, index=self.data.index, columns=self.data.columns)
                count_matrix = pd.DataFrame(0, index=self.data.index, columns=self.data.columns)

                for matrix in matrices:
                    count_matrix += ~matrix.data.isna()
                    combo_matrix += matrix.data.fillna(0)

                self.data = combo_matrix / count_matrix

            #def print_matrix(self):
                #print(self.data)

        ##### create Matrix objects from the expanded_matrices list
        partial_matrices = [Matrix(matrix) for matrix in expanded_matrices]

        ##### display partial matrices:
        #for index, partial_matrix in enumerate(partial_matrices):
            #print(f"Partial Matrix {index + 1}:")
            #partial_matrix.print_matrix()
            #print()

        ##### generate an object-oriented scratch matrix from partial matrices
        class ScratchMatrix:
            def __init__(self, rows, columns):
                self.rows = rows
                self.columns = columns
                self.matrix = pd.DataFrame(0, index=rows, columns=columns)
                self.count_matrix = pd.DataFrame(0, index=rows, columns=columns)

            def add_partial_matrix(self, partial_matrix):
                self.count_matrix += ~partial_matrix.data.isna()
                self.matrix += partial_matrix.data.fillna(0)    

            def compute_average_matrix(self):
                self.matrix = self.matrix / self.count_matrix

        #def display_matrix(self):
            #print("scratch_matrix:")
            #print(self.matrix)

        ##### initialize the scratch matrix object
        scratch_matrix = ScratchMatrix(tip_labels_list, tip_labels_list)

        ####################################################################################################
        ##### loop through partial matrices ################################################################
        ####################################################################################################

        ##### loop through the partial matrices
        for partial_matrix in partial_matrices:
            scratch_matrix.add_partial_matrix(partial_matrix)

        ##### compute the cumulative average matrix
        scratch_matrix.compute_average_matrix()

        ##### display the scratch matrix
        #scratch_matrix.display_matrix()

        ##### find the smallest distance and the sister pair in the scratch matrix
        min_distance = np.nanmin(scratch_matrix.matrix.values)
        min_indices = np.where(scratch_matrix.matrix.values == min_distance)
        sister_pair = (scratch_matrix.matrix.index[min_indices[0][0]], scratch_matrix.matrix.columns[min_indices[1][0]])

        ##### create the composite cluster name
        composite_cluster = f"{sister_pair[0]},{sister_pair[1]}"

        ##### compute distances between the composite cluster and the rest of the leaf nodes
        composite_distances = {}
        for leaf in scratch_matrix.matrix.index:
            if leaf != sister_pair[0] and leaf != sister_pair[1]:
                composite_distances[leaf] = (scratch_matrix.matrix.loc[sister_pair[0], leaf] + scratch_matrix.matrix.loc[sister_pair[1], leaf]) / 2

        ##### remove the sister pair from the scratch matrix
        condensed_matrix = scratch_matrix.matrix.drop(index=[sister_pair[0], sister_pair[1]], columns=[sister_pair[0], sister_pair[1]])

        ##### add the composite cluster to the condensed matrix
        condensed_matrix = pd.concat([condensed_matrix, pd.Series(composite_distances, name=composite_cluster)], axis=1)
        condensed_matrix.loc[composite_cluster] = pd.Series(composite_distances)
        condensed_matrix.loc[composite_cluster, composite_cluster] = np.nan

        ##### find the indices of the identified cluster in each partial matrix
        for partial_matrix in partial_matrices:
            if sister_pair[0] in partial_matrix.data.index and sister_pair[1] in partial_matrix.data.columns:
                row_index = partial_matrix.data.index.get_loc(sister_pair[0])
                col_index = partial_matrix.data.columns.get_loc(sister_pair[1])
                partial_matrix.data.iat[row_index, col_index] = np.nan
                partial_matrix.data.iat[col_index, row_index] = np.nan

        ##### add the composite cluster and its distances to each partial matrix
        for partial_matrix in partial_matrices:
            if composite_cluster not in partial_matrix.data.index:
                partial_matrix.data = pd.concat([partial_matrix.data, pd.Series(np.nan, name=composite_cluster, index=partial_matrix.data.index)], axis=1)
            if composite_cluster not in partial_matrix.data.columns:
                partial_matrix.data.loc[composite_cluster] = pd.Series(np.nan, index=partial_matrix.data.columns)

            for leaf, distance in composite_distances.items():
                partial_matrix.data.loc[composite_cluster, leaf] = distance
                partial_matrix.data.loc[leaf, composite_cluster] = distance

        ##### remove the clustered elements from the partial_matrices
        for partial_matrix in partial_matrices:
            partial_matrix.data = partial_matrix.data.drop([sister_pair[0], sister_pair[1]], axis=0)
            partial_matrix.data = partial_matrix.data.drop([sister_pair[0], sister_pair[1]], axis=1)

        ##### create a dataframe to store cluster values
        clusters = pd.DataFrame(columns=["lineage_a", "lineage_b", "pairwise_distance"])

        ##### find the sister pair and pairwise distance
        sister_pair_distance = scratch_matrix.matrix.loc[sister_pair[0], sister_pair[1]]

        ##### update the clusters DataFrame
        clusters.loc[len(clusters)] = [sister_pair[0], sister_pair[1], sister_pair_distance]

        while len(condensed_matrix) > 1:
            ##### find the smallest distance and the sister pair in the condensed matrix
            min_distance = np.nanmin(condensed_matrix.values)
            min_indices = np.where(condensed_matrix.values == min_distance)
            #print("condensed_matrix:")
            #print(condensed_matrix)
            sister_pair = (condensed_matrix.index[min_indices[0][0]], condensed_matrix.columns[min_indices[1][0]])

            ##### create the composite cluster name
            composite_cluster = f"{sister_pair[0]},{sister_pair[1]}"

            ##### compute distances between the composite cluster and the rest of the leaf nodes
            composite_distances = {}
            for leaf in condensed_matrix.index:
                if leaf != sister_pair[0] and leaf != sister_pair[1]:
                    if sister_pair[0] in condensed_matrix.index and leaf in condensed_matrix.columns:
                        distance_a = condensed_matrix.loc[sister_pair[0], leaf]
                    else:
                        distance_a = np.nan

                    if sister_pair[1] in condensed_matrix.index and leaf in condensed_matrix.columns:
                        distance_b = condensed_matrix.loc[sister_pair[1], leaf]
                    else:
                        distance_b = np.nan

                    if not np.isnan(distance_a) and not np.isnan(distance_b):
                        composite_distances[leaf] = (distance_a + distance_b) / 2
                    elif not np.isnan(distance_a):
                        composite_distances[leaf] = distance_a
                    elif not np.isnan(distance_b):
                        composite_distances[leaf] = distance_b

            ##### remove the sister pair from the condensed matrix
            condensed_matrix = condensed_matrix.drop(index=[sister_pair[0], sister_pair[1]], columns=[sister_pair[0], sister_pair[1]])

            ##### add the composite cluster to the condensed matrix
            condensed_matrix = pd.concat([condensed_matrix, pd.Series(composite_distances, name=composite_cluster)], axis=1)
            condensed_matrix.loc[composite_cluster] = pd.Series(composite_distances)
            condensed_matrix.loc[composite_cluster, composite_cluster] = np.nan

            ##### find the indices of the identified cluster in each partial matrix
            for partial_matrix in partial_matrices:
                if sister_pair[0] in partial_matrix.data.index and sister_pair[1] in partial_matrix.data.columns:
                    row_index = partial_matrix.data.index.get_loc(sister_pair[0])
                    col_index = partial_matrix.data.columns.get_loc(sister_pair[1])
                    partial_matrix.data.iat[row_index, col_index] = np.nan
                    partial_matrix.data.iat[col_index, row_index] = np.nan

            ##### add the composite cluster and its distances to each partial matrix
            for partial_matrix in partial_matrices:
                if composite_cluster not in partial_matrix.data.index:
                    partial_matrix.data = pd.concat([partial_matrix.data, pd.Series(np.nan, name=composite_cluster, index=partial_matrix.data.index)], axis=1)
                if composite_cluster not in partial_matrix.data.columns:
                    partial_matrix.data.loc[composite_cluster] = pd.Series(np.nan, index=partial_matrix.data.columns)

                for leaf, distance in composite_distances.items():
                    partial_matrix.data.loc[composite_cluster, leaf] = distance
                    partial_matrix.data.loc[leaf, composite_cluster] = distance

            ##### remove the clustered elements from the partial_matrices
            for partial_matrix in partial_matrices:
                partial_matrix.data = partial_matrix.data.drop([sister_pair[0], sister_pair[1]], axis=0)
                partial_matrix.data = partial_matrix.data.drop([sister_pair[0], sister_pair[1]], axis=1)

            ##### store cluster information in the clusters DataFrame
            clusters.loc[len(clusters)] = [sister_pair[0], sister_pair[1], min_distance]

        ##### Export clusters as a CSV file
        clusters.to_csv("clusters_and_pairwise_distances_list.csv", index=False)

        ##### Display the list of clusters and pairwise distances
        #print("clusters_and_paiwise_distances:")
        #print(clusters)

        ##### Convert clusters into a square pairwise distance matrix
        df = pd.DataFrame(clusters, columns=['lineage_a', 'lineage_b', 'pairwise_distance'])

        lineages = set(df['lineage_a'].str.split(',').sum() + df['lineage_b'].str.split(',').sum())
        lineages = sorted(lineages)

        final_pairwise_matrix = pd.DataFrame(np.nan, index=lineages, columns=lineages)

        for _, row in df.iterrows():
            lineage_a = row['lineage_a']
            lineage_b = row['lineage_b']
            distance = row['pairwise_distance']

            lineages_a = lineage_a.split(',')
            lineages_b = lineage_b.split(',')

            for a in lineages_a:
                for b in lineages_b:
                    final_pairwise_matrix.at[a, b] = distance
                    final_pairwise_matrix.at[b, a] = distance

        final_pairwise_matrix = final_pairwise_matrix.astype(float)

        ##### Display the list of clusters and pairwise distances
        #print("final_pairwise_matrix:")
        #print(final_pairwise_matrix)

        ##### Export final pairwise matrix as a CSV file
        final_pairwise_matrix.to_csv("final_pairwise_distance_matrix.csv", index=True, header=True, na_rep='NaN')

        ####################################################################################################
        ##### part04 #######################################################################################
        ##### construct UPGMA supertree from final pairwise matrix ########################################
        ####################################################################################################

        ##### Set diagonal of final pairwise matrix to 0
        final_pairwise_matrix_no_nans = final_pairwise_matrix.copy()
        np.fill_diagonal(final_pairwise_matrix_no_nans.values, 0)

        ##### Convert the upper triangular matrix to a condensed distance matrix
        distances = squareform(final_pairwise_matrix_no_nans)
        distances = distances / 2

        ##### Perform UPGMA clustering
        upgma = hierarchy.linkage(distances, method='average')

        ##### Plot the dendrogram
        plt.figure(figsize=(10, 8))
        dendrogram = hierarchy.dendrogram(upgma, labels=final_pairwise_matrix_no_nans.columns.tolist())
        plt.xlabel('Samples', fontsize=12)
        plt.ylabel('Time', fontsize=12)
        plt.title('Chrono-STA supertimetree', fontsize=14)
        plt.tight_layout()
        #plt.show()

        ##### Define function
        def get_newick(node, parent_dist, leaf_names, newick=''):
            """
            Convert scipy.cluster.hierarchy.to_tree()-output to Newick format.
            :param node: output of scipy.cluster.hierarchy.to_tree()
            :param parent_dist: distance of the parent node
            :param leaf_names: list of leaf names
            :param newick: leave empty, this variable is used in recursion.
            :returns: tree in Newick format
            """
            if node.is_leaf():
                return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
            else:
                if len(newick) > 0:
                    newick = "):%.2f%s" % (parent_dist - node.dist, newick)
                else:
                    newick = ");"
                newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
                newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
                newick = "(%s" % (newick)
                return newick

        ##### Define leaf names from matrix columns
        leaf_names = final_pairwise_matrix_no_nans.columns.tolist()

        ##### Convert the linkage matrix to a tree object
        tree = hierarchy.to_tree(upgma)

        ##### Convert the tree to Newick format
        newick_tree = get_newick(tree, tree.dist, leaf_names)
        #print("UPGMA supertimetree from pairwise distance matrix")
        #print(newick_tree)

        ##### Export the Newick tree as a file
        with open('supertimetree_from_final_pairwise_distance_matrix.newick', 'w') as file:
            file.write(newick_tree)

        ##### Load the tree from a file (if required)
        #newick_tree = Phylo.read("supertimetree_from_final_pairwise_distance_matrix.newick", "newick")

        ##### Draw phylogenetic tree
        #Phylo.draw(newick_tree, do_show=True)

        ####################################################################################################
        ##### make supertimetree ultrametric ###############################################################
        ####################################################################################################

        ##### Import necessary R packages
        phytools = importr("phytools")

        ##### Read the tree from the Newick file
        tree = robjects.r['read.tree']("supertimetree_from_final_pairwise_distance_matrix.newick")

        ##### Ladderize the tree
        ape = importr("ape")
        ladderized_tree = ape.ladderize(tree)

        ##### Force the tree to be ultrametric using NNLS method
        phy2 = phytools.force_ultrametric(ladderized_tree, method="nnls", message="FALSE")

        ##### Write the ultrametric tree to a new Newick file
        robjects.r['write.tree'](phy2, file="chronosta_supertimetree.newick")

        ##### Read the Newick representation of the ultrametric tree
        with open("chronosta_supertimetree.newick", "r") as f:
            newick_tree2 = f.read()

        print("Chrono-STA supertimetree:")
        print(newick_tree2)

        ##### Generate pairwise distance matrix from chronosta_supertimetree.newick
        def parse_newick_tree_from_file(newick_file):
            """
            Parse Newick tree from file
            """
            with open(newick_file, 'r') as f:
                newick_str = f.readline().strip()
            return Tree(newick_str)

        def compute_pairwise_distances(tree):
            """
            Compute pairwise distances between clusters
            """
            pairwise_distances = []
            for leaf1 in tree:
                for leaf2 in tree:
                    if leaf1 != leaf2:
                        distance = tree.get_distance(leaf1, leaf2)
                        pairwise_distances.append((leaf1.name, leaf2.name, distance))
            return pairwise_distances

        def format_table(pairwise_distances):
            """
            Format pairwise distances into a table
            """
            table = "lineage_a,lineage_b,p            airwise_distance\n"
            for pair in pairwise_distances:
                table += f"{pair[0]},{pair[1]},{pair[2]}\n"
            return table

        def clusters_to_pairwise_matrix(clusters):
            lineages = set()
            for index, row in clusters.iterrows():
                lineages.update(row['lineage_a'].split(','))
                lineages.update(row['lineage_b'].split(','))

            lineages = sorted(lineages)
            chronosta_pairwise_matrix = pd.DataFrame(np.nan, index=lineages, columns=lineages)

            for index, row in clusters.iterrows():
                lineage_a = row['lineage_a']
                lineage_b = row['lineage_b']
                distance = row['pairwise_distance']

                lineages_a = lineage_a.split(',')
                lineages_b = lineage_b.split(',')

                for a in lineages_a:
                    for b in lineages_b:
                        chronosta_pairwise_matrix.at[a, b] = distance
                        chronosta_pairwise_matrix.at[b, a] = distance

            return chronosta_pairwise_matrix

        if __name__ == "__main__":

            ##### Parse Newick tree from file
            newick_file = "chronosta_supertimetree.newick"
            tree = parse_newick_tree_from_file(newick_file)

            ##### Compute pairwise distances between clusters
            pairwise_distances = compute_pairwise_distances(tree)

            ##### Format pairwise distances into a table
            table = format_table(pairwise_distances)
            #print("Pairwise Distance Table:")
            #print(table)

            ##### Convert clusters into a square pairwise distance matrix
            df = pd.DataFrame(pairwise_distances, columns=['lineage_a', 'lineage_b', 'pairwise_distance'])
            chronosta_pairwise_matrix = clusters_to_pairwise_matrix(df)

            ##### Display the pairwise distance matrix
            #print("\Chrono-STA Pairwise Distance Matrix:")
            #print(chronosta_pairwise_matrix)

            ##### Export final pairwise distance matrix as a CSV file
            chronosta_pairwise_matrix.to_csv("chronosta_supertimetree_pairwise_distance_matrix.csv", index=True, header=True, na_rep='NaN')
            
    except Exception as e:
        warnings.warn(f"{str(e)}")
        print("ERROR: An error has ocurred. Ensure all subtrees or subsets within the set have common tip labels with the full set.")

if __name__ == "__main__":
    main()

####################################################################################################
##### end of the algorithm #########################################################################
####################################################################################################
