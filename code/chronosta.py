#!/usr/bin/python3

###################################################################################
##### Chronological Supertree Algorithm (Chrono-STA) ##############################
##### Jose Barba, Jack Craig, and Sudhir Kumar ####################################
###################################################################################

# ────────────────────────────────────────────────────────────────────────────────
# print heading
# ────────────────────────────────────────────────────────────────────────────────

print("\n\n***** Chronological Supertree Algorithm (Chrono-STA 1.3 built April 25 2025)\n"
      "***** Developed by Jose Barba, Jack Craig and Sudhir Kumar\n")

# ────────────────────────────────────────────────────────────────────────────────
# set up python env
# ────────────────────────────────────────────────────────────────────────────────

import os
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────────────────────────────────────
# helper classes
# ────────────────────────────────────────────────────────────────────────────────

class Matrix:
    def __init__(self, data: pd.DataFrame):
        self.data = data

class ScratchMatrix:
    def __init__(self, labels: List[str]):
        self.matrix = pd.DataFrame(0.0, index=labels, columns=labels)
        self.count = pd.DataFrame(0,   index=labels, columns=labels)

    def add(self, partial: Matrix):
        self.count  += ~partial.data.isna()
        self.matrix += partial.data.fillna(0.0)

    def mean(self):
        self.matrix = self.matrix / self.count

# ────────────────────────────────────────────────────────────────────────────────
# build per‑tree pairwise distance matrices
# ────────────────────────────────────────────────────────────────────────────────

def read_trees() -> Tuple[List[pd.DataFrame], List[str]]:
    files = [f for f in os.listdir('.') if f.endswith('.nwk')]
    tips: set[str] = set()
    for f in files:
        for cl in Phylo.read(f, 'newick').get_terminals():
            tips.add(cl.name)
    tips = sorted(tips)

    mats: List[pd.DataFrame] = []
    for f in files:
        t = Phylo.read(f, 'newick')
        pairs = {(a.name, b.name): t.distance(a, b)
                 for a in t.get_terminals() for b in t.get_terminals() if a != b}
        df = pd.DataFrame.from_dict(pairs, orient='index', columns=['d'])
        df[['r', 'c']] = pd.DataFrame(df.index.tolist(), index=df.index)
        mat = df.pivot(index='r', columns='c', values='d')
        mats.append(mat)
    print("\nInput data:")
    print(f"Number of trees in the set: {len(files)}")
    print(f"Number of unique tip labels: {len(tips)}\n")
    return mats, tips

# ────────────────────────────────────────────────────────────────────────────────


def expand(mats: List[pd.DataFrame], labels: List[str]) -> List[pd.DataFrame]:
    out = []
    for m in mats:
        e = pd.DataFrame(np.nan, index=labels, columns=labels)
        common_r = m.index.intersection(labels)
        common_c = m.columns.intersection(labels)
        e.loc[common_r, common_c] = m.loc[common_r, common_c]
        out.append(e)
    return out

# ────────────────────────────────────────────────────────────────────────────────
# weighted update helper
# ────────────────────────────────────────────────────────────────────────────────

def wavg(dik: float, djk: float, ni: int, nj: int) -> float:
    if np.isnan(dik) and np.isnan(djk):
        return np.nan
    if np.isnan(dik):
        return djk
    if np.isnan(djk):
        return dik
    return (ni * dik + nj * djk) / (ni + nj)

# ────────────────────────────────────────────────────────────────────────────────
# parallel weighted clustering loop
# ────────────────────────────────────────────────────────────────────────────────

def build_supermatrix(expanded: List[pd.DataFrame], labels: List[str]) -> pd.DataFrame:
    partials = [Matrix(m.copy()) for m in expanded]

    # cluster sizes
    csize: Dict[str, int] = {l: 1 for l in labels}

    scratch = ScratchMatrix(labels)
    for p in partials:
        scratch.add(p)
    scratch.mean()

    condensed = scratch.matrix.copy()
    clusters = []

    while len(condensed) > 1:
        md = np.nanmin(condensed.values)
        i_idx, j_idx = np.where(condensed.values == md)
        i = condensed.index[i_idx[0]]
        j = condensed.columns[j_idx[0]]

        clusters.append({'lineage_a': i, 'lineage_b': j, 'pairwise_distance': md})

        comp = f"{i},{j}"
        ni, nj = csize[i], csize[j]
        comp_d: Dict[str, float] = {}
        for k in condensed.index:
            if k in (i, j):
                continue
            dik = condensed.at[i, k] if (i in condensed.index and k in condensed.columns) else np.nan
            djk = condensed.at[j, k] if (j in condensed.index and k in condensed.columns) else np.nan
            comp_d[k] = wavg(dik, djk, ni, nj)

        # remove old rows/cols
        condensed = condensed.drop(index=[i, j], columns=[i, j])
        # add new col first (avoids empty‑frame error)
        condensed[comp] = pd.Series(comp_d)
        condensed.loc[comp] = pd.Series(comp_d)
        condensed.at[comp, comp] = np.nan

        # update cluster size dict
        csize[comp] = ni + nj
        del csize[i], csize[j]

        # propagate into every partial matrix
        for p in partials:
            # set i‑j to NaN where present
            if i in p.data.index and j in p.data.columns:
                p.data.at[i, j] = p.data.at[j, i] = np.nan
            # drop i and j
            p.data = p.data.drop(index=[x for x in (i, j) if x in p.data.index], errors='ignore')
            p.data = p.data.drop(columns=[x for x in (i, j) if x in p.data.columns], errors='ignore')
            # ensure comp present (col first, then row)
            if comp not in p.data.columns:
                p.data[comp] = np.nan
            if comp not in p.data.index:
                p.data.loc[comp] = np.nan
            # fill distances
            for k, d in comp_d.items():
                if k in p.data.columns:
                    p.data.at[comp, k] = p.data.at[k, comp] = d

    return pd.DataFrame(clusters)

# ────────────────────────────────────────────────────────────────────────────────
# convert clusters into pairwise matrix
# ────────────────────────────────────────────────────────────────────────────────

def clusters_to_matrix(cl: pd.DataFrame) -> pd.DataFrame:
    taxa: set[str] = set()
    for _, r in cl.iterrows():
        taxa.update(r['lineage_a'].split(','))
        taxa.update(r['lineage_b'].split(','))
    taxa = sorted(taxa)
    mat = pd.DataFrame(np.nan, index=taxa, columns=taxa)
    for _, r in cl.iterrows():
        a_l = r['lineage_a'].split(',')
        b_l = r['lineage_b'].split(',')
        d = r['pairwise_distance']
        for a in a_l:
            for b in b_l:
                mat.at[a, b] = mat.at[b, a] = d
    return mat

# ────────────────────────────────────────────────────────────────────────────────
# build and save final supertree
# ────────────────────────────────────────────────────────────────────────────────

def build_supertree(mat: pd.DataFrame, prefix: str = 'chronosta'):
    m = mat.copy()
    np.fill_diagonal(m.values, 0.0)
    # Convert patristic tip‑to‑tip distances to divergence time (halve)
    dist_vec = squareform(m) / 2.0  # patristic → divergence
    link = hierarchy.linkage(dist_vec, method='average' )

    plt.figure(figsize=(8, 6))
    hierarchy.dendrogram(link, labels=m.index.tolist())
    plt.tight_layout()
    plt.savefig(f'{prefix}_dendrogram.png')

    def to_newick(node, parent, names, nw=''):
        if node.is_leaf():
            return f"{names[node.id]}:{parent - node.dist:.6f}{nw}"
        if nw:
            nw = f"):{parent - node.dist:.6f}{nw}"
        else:
            nw = ');'
        nw = to_newick(node.get_left(), node.dist, names, nw)
        nw = to_newick(node.get_right(), node.dist, names, f",{nw}")
        return f"({nw}"

    tree = hierarchy.to_tree(link)
    newick = to_newick(tree, tree.dist, m.index.tolist())
    with open(f'{prefix}_supertree.newick', 'w') as fh:
        fh.write(newick)

    # Read and print the Newick tree to display it on the screen
    with open(f'{prefix}_supertree.newick', 'r') as fh:
        newick_content = fh.read()
    print("\nChrono‑STA supertree:")
    print(newick_content)

# ────────────────────────────────────────────────────────────────────────────────
# main
# ────────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    partial_mats, tips = read_trees()
    expanded_mats = expand(partial_mats, tips)
    clusters_df = build_supermatrix(expanded_mats, tips)
    clusters_df.to_csv('clusters_and_pairwise_distances_list.csv', index=False)

    final_mat = clusters_to_matrix(clusters_df)
    final_mat.to_csv('final_pairwise_distance_matrix.csv')

    build_supertree(final_mat)

# ────────────────────────────────────────────────────────────────────────────────
# final notes
# ────────────────────────────────────────────────────────────────────────────────

    print('\n\nResults written to:')
    print('List of clusters and pairwise distances: ' 'clusters_and_pairwise_distances_list.csv')
    print('Final supermatrix: ' 'final_pairwise_distance_matrix.csv')
    print('Chrono‑STA supertree in newick format: ' 'chronosta_supertree.newick')
    print('Chrono‑STA supertree in png format: ' 'chronosta_dendrogram.png\n\n')
