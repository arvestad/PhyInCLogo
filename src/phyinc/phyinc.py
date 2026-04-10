import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import logomaker
import pandas as pd

from Bio import Phylo
from io import StringIO

import phyinc.io as io
from phyinc.colorhelper import random_color_scheme


_seq_type_map = {
    'dna':     io.unambiguous_dna_alphabet,
    'rna':     io.unambiguous_rna_alphabet,
    'protein': io.unambiguous_protein_alphabet,
}

_default_color_scheme = {
    'dna':     'classic',
    'rna':     'classic',
    'protein': 'chemistry',
}


def create_logo(alignment, tree, seq_type, color_scheme='guess', bar_width=None, max_bar_height=None, title=None):
    """
    Compute a PIC-weighted sequence logo and return it as a matplotlib Figure.

    Parameters
    ----------
    alignment : Bio.Align.MultipleSeqAlignment
        Multiple sequence alignment; record IDs must match tree leaf labels.
    tree : Bio.Phylo.BaseTree.Tree
        Rooted phylogenetic tree.
    seq_type : str
        One of 'dna', 'rna', or 'protein'.
    color_scheme : str
        Logomaker color scheme string, or 'guess' to pick based on seq_type.
        Default: 'guess'.
    bar_width : float, optional
        Width of each position column in inches.
    max_bar_height : float, optional
        Y-axis ceiling in bits.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if seq_type not in _seq_type_map:
        raise ValueError(f"seq_type must be one of {list(_seq_type_map)}, got {seq_type!r}")
    alphabet = _seq_type_map[seq_type]

    if color_scheme == 'guess':
        color_scheme = _default_color_scheme[seq_type]
    elif color_scheme == 'random':
        color_scheme = random_color_scheme(alphabet.letters())

    if isinstance(alignment, dict):
        seq_dict = alignment
        seq_length = len(next(iter(seq_dict.values())).seq)
    else:
        seq_dict = {record.id: record for record in alignment}
        seq_length = alignment.get_alignment_length()

    array = compute_pic_array(tree, seq_dict, alphabet, seq_length)

    df = pd.DataFrame(array, columns=list(alphabet.letters()))

    # Schneider & Stephens (1990) gap correction:
    # row sums give the PIC-weighted non-gap fraction at each position.
    non_gap_weight = df.sum(axis=1)
    df = df.div(non_gap_weight, axis=0)          # conditional probabilities (given non-gap)
    df = logomaker.transform_matrix(df, from_type='probability', to_type='information')
    df = df.mul(non_gap_weight, axis=0)          # scale IC by non-gap fraction

    col_width = bar_width if bar_width is not None else 0.4
    fig_width = max(4, seq_length * col_width)
    fig, ax = plt.subplots(figsize=(fig_width, 2.5))

    logomaker.Logo(df, ax=ax, color_scheme=color_scheme)
    ax.set_xlabel('Position')
    ax.set_ylabel('Information content (bits)')

    if max_bar_height is not None:
        ax.set_ylim(0, max_bar_height)

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    return fig


def convert_matrix_to_array(matrix, seq_type, seq_length):
    return matrix.T



def collapse_bifurcations(children, seq_dict, zero_matrix):
    """
    TODO: misleading name! I completely misunderstood Haolin's code.

    In case the tree has multifurcating nodes, we create bifurcations
    by adding very short edges.
    """
    branch_length_temp = []
    seq_matrix_temp = []
    for child in children:
        branch_length_temp.append(float(child.branch_length))
        seq_matrix_temp.append(seq_dict[child])

    while len(branch_length_temp) > 1 and len(seq_matrix_temp) > 1:
        leaf_i_branch_length = branch_length_temp.pop(0)
        leaf_j_branch_length = branch_length_temp.pop(0)
        leaf_i_seq_matrix = seq_matrix_temp.pop(0)
        leaf_j_seq_matrix = seq_matrix_temp.pop(0)

        if (leaf_i_branch_length == 0) and (leaf_j_branch_length == 0):
            branch_length_temp.append(0)
            seq_matrix_temp.append(zero_matrix)
        else:
            branch_length_temp.append(
                (
                    (leaf_i_branch_length * leaf_j_branch_length)
                    / (leaf_i_branch_length + leaf_j_branch_length)
                )
            )
            seq_matrix_temp.append(
                np.add(
                    (
                        leaf_i_branch_length
                        / (leaf_i_branch_length + leaf_j_branch_length)
                    )
                    * leaf_i_seq_matrix,
                    (
                        leaf_j_branch_length
                        / (leaf_i_branch_length + leaf_j_branch_length)
                    )
                    * leaf_j_seq_matrix,
                )
            )
    return branch_length_temp[0], seq_matrix_temp[0]


def compute_pic_array(tree, alignment, seq_type, seq_length):
    """
    Run the PIC traversal and return the frequency array with shape
    (seq_length, n_chars), ready for use with Logomaker or export.
    """
    zero_matrix = np.zeros((len(seq_type), seq_length), dtype=float)
    seq_dict = {}
    for child in tree.clade:
        traverse_postorder(child, alignment, seq_type, zero_matrix, seq_dict)
    _, seq_matrix = collapse_bifurcations(tree.clade, seq_dict, zero_matrix)
    return convert_matrix_to_array(seq_matrix, seq_type, seq_length)


def traverse_postorder(clade, alignment, seq_type, zero_matrix, seq_dict):
    if len(clade) == 0:  # leaf
        seq = str(alignment[clade.name].seq)
        seq_matrix = zero_matrix.copy()
        for i, char in enumerate(seq.upper()):
            if char not in ('-', '.'):
                seq_matrix[seq_type.letters().index(char), i] = 1.0
        seq_dict[clade] = seq_matrix
    else:  # internal node
        for child in clade:
            traverse_postorder(child, alignment, seq_type, zero_matrix, seq_dict)
        branch_length, seq_matrix = collapse_bifurcations(clade, seq_dict, zero_matrix)
        clade.branch_length = float(clade.branch_length) + branch_length
        seq_dict[clade] = seq_matrix
