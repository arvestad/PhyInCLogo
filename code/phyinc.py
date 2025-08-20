import sys
import os
import time
import math
import argparse
import warnings
import config

import numpy as np
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from pathlib import Path
from Bio import Phylo
from Bio import SeqIO


from weblogo import LogoData, LogoOptions, LogoFormat, png_formatter, pdf_formatter
from weblogo.seq import (
    generic_alphabet,
    protein_alphabet,
    nucleic_alphabet,
    dna_alphabet,
    rna_alphabet,
    reduced_nucleic_alphabet,
    reduced_protein_alphabet,
    unambiguous_dna_alphabet,
    unambiguous_rna_alphabet,
    unambiguous_protein_alphabet,
)


def export_scores_to_file(scores_str, file_name):
    try:
        with open(file_name, "x") as f:
            f.write(scores_str)
    except FileExistsError:
        print("Already exists.")


"""
def convert_matrix_to_scores(matrix):
    #convert the seq_counter matrix to the format used by pyseqlogo 
    for c in config.charaters:
        if sum(config.seq_counter[c]) == 0 or sum(config.seq_counter[c]) == float(0):
            config.seq_counter.pop(c)
        elif c == 'X': # removes counts for X in sequnces since it indicates any/unknown
            config.seq_counter.pop(c)
        else:
            config.existing_characters.append(c)

    all_scores = []
    for i in range(0,config.seq_length):
        scores = []
        for c in config.existing_characters:
            scores.append((str(c),matrix[config.charaters.index(c)][i]))
        all_scores.append(scores)
    export_scores_to_file(str(all_scores),'ex1_t3.txt')
    return all_scores

def convert_count_to_scores(matrix):
    #convert the count dictionary of original data (sequnces) to the format used by pyseqlogo 
    all_scores_count = []
    for i in range(0,config.seq_length):
        position_scores = []
        for key, value in matrix.items():
            position_scores.append((str(key),float(value[i]/config.terminals)))
        all_scores_count.append(position_scores)
    export_scores_to_file(str(all_scores_count),'ex1_without.txt')
    return all_scores_count

"""


def convert_matrix_to_array(matrix):
    array = []
    for i in range(0, config.seq_length):
        sub_array = []
        for c in config.seq_type:
            sub_array.append(matrix[config.seq_type.letters().index(c)][i])
        array.append(sub_array)
    # export_scores_to_file(str(all_scores),'ex1_t3.txt')
    return np.array(array)


def convert_count_to_array(matrix):
    array = []
    for i in range(0, config.seq_length):
        sub_array = []
        for c in config.seq_type:
            if str(c) not in matrix:
                sub_array.append(0)
            else:
                sub_array.append(matrix[c][i])
        array.append(sub_array)
    # export_scores_to_file(str(all_scores_count),'ex1_without.txt')
    return np.array(array)


def set_seq(leaf_i, leaf_j):
    """caculate frequncy matrix (x_i) for parent node"""
    l = float(leaf_i.branch_length)
    r = float(leaf_j.branch_length)
    if (l == 0) and (r == 0):
        return config.matrix
    else:
        seq_matrix = np.add(
            (l / (l + r)) * config.seq_dict[leaf_j],
            (r / (l + r)) * config.seq_dict[leaf_i],
        )
        return seq_matrix


def add_length(leaf_i, leaf_j):
    """caculate new branch length (v_i') for parent node"""
    l = float(leaf_i.branch_length)
    r = float(leaf_j.branch_length)
    if (l == 0) and (r == 0):
        return float(0)
    else:
        return (l * r) / (l + r)


def unbifrucating(childs):
    """For case where a node has more than 2 childs"""
    branch_length_temp = []
    seq_matrix_temp = []
    for child in childs:
        branch_length_temp.append(float(child.branch_length))
        seq_matrix_temp.append(config.seq_dict[child])

    while len(branch_length_temp) > 1 and len(seq_matrix_temp) > 1:
        leaf_i_branch_length = branch_length_temp.pop(0)
        leaf_j_branch_length = branch_length_temp.pop(0)
        leaf_i_seq_matrix = seq_matrix_temp.pop(0)
        leaf_j_seq_matrix = seq_matrix_temp.pop(0)

        if (leaf_i_branch_length == 0) and (leaf_j_branch_length == 0):
            branch_length_temp.append(0)
            seq_matrix_temp.append(config.matrix)
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


def PIC_postorder(tree):

    for child in tree.clade:
        traverse_postorder(child)

    result, seq_matrix = unbifrucating(tree.clade)

    array = convert_matrix_to_array(seq_matrix)
    print(array)

    logo_data = LogoData.from_counts(alphabet=config.seq_type, counts=array)

    logo_options = LogoOptions()
    logo_options.title = "With PIC logo"
    logo_options.stack_width = 50  # increase width of each position
    logo_options.stack_height = 100  # increase overall height

    logo_format = LogoFormat(logo_data, logo_options)

    # Save as PNG
    with open("With_PIC_logo.png", "wb") as f:
        f.write(png_formatter(logo_data, logo_format))

    """
    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(ALL_SCORES1,seq_type='aa', yaxis='probability', colorscheme='hydrophobicity')
    fig.tight_layout()
    fig.show()
    fig.savefig("with_PIC.png")
    print("with_PIC.png saved")
    

    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(ALL_SCORES1,data_type='bits',seq_type=config.seq_type, yaxis='bits',colorscheme='hydrophobicity',draw_axis=True)
    fig.tight_layout()
    fig.show()
    fig.savefig("with_PIC_bits.png")
    print("with_PIC_bits.png saved")
    """


def traverse_postorder(clade):
    if len(clade) == 0:  # only tips of the tree will have length 0
        clade.seq = str(
            config.updated_dict[clade.name].seq
        )  # store str(sequnces) as an artribute for the clade object(of biopython package).

        seq_matrix = config.matrix.copy()  # make a copy of the default sequnce matrix
        for i in range(0, len(clade.seq)):
            config.seq_counter[clade.seq[i].upper()][
                i
            ] += 1  # count and store for each character as a whole

            character_index = config.seq_type.letters().index(clade.seq[i].upper())
            seq_matrix[character_index, i] = float(
                1
            )  # a count matrix for each individual leaf, used later to caculate frequncy matrix for parent nodes

        config.seq_dict[clade] = seq_matrix  # stores the matrix to dictionary

    if len(clade) > 0:  # a parent node
        for child in clade:
            traverse_postorder(child)
        if len(clade) == 2:
            clade.branch_length = float(clade.branch_length) + add_length(
                clade[0], clade[1]
            )
            config.seq_dict[clade] = set_seq(clade[0], clade[1])
        if len(clade) > 2:
            branch_length, seq_matrix = unbifrucating(clade)
            clade.branch_length = float(clade.branch_length) + branch_length
            config.seq_dict[clade] = seq_matrix


def parse_dict(filename, filetype="fasta"):
    """checks length and type of sequnces."""

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    updated_dict = {}

    for key, value in record_dict.items():
        if len(value) != config.seq_length:
            if config.seq_length == 0:
                config.seq_length = len(value)
            else:
                raise Exception(
                    "Sequence length from provided fastafile are inconsistent"
                )

        config.charaters.update(set(value.upper()))

        updated_dict[key] = value

    return updated_dict


def find_clades(clade, condition):
    """Find clades matching a condition, used in testing."""
    matches = []
    if condition(clade):
        matches.append(clade)
    for sub_clade in clade.clades:
        matches.extend(find_clades(sub_clade, condition))
    return matches


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Apply PIC on tree file and corresbond Name-sequnce file, outputs 2 .png file(sequnce logo with and with out applying PIC)."
    )
    parser.add_argument(
        "input_tree", type=str, help="Path to the .tree file (newick format)"
    )
    parser.add_argument(
        "input_fa",
        type=str,
        help="Path to the .fa file (fasta format), assumes sequnce name in the format of: ETA_STAAU/96-110",
    )
    args = parser.parse_args()

    # Validate file paths
    tree_file = Path(args.input_tree)
    fa_file = Path(args.input_fa)

    if not tree_file.is_file():
        raise FileNotFoundError(f".tree {tree_file} does not exist.")
    if not fa_file.is_file():
        raise FileNotFoundError(f".fa {fa_file} does not exist.")

    print(f"Processing: {tree_file} and {fa_file}")

    tree = Phylo.read(tree_file, "newick")
    config.terminals = tree.count_terminals()

    # config.py to store global variables
    config.updated_dict = parse_dict(fa_file, "fasta")
    print("sequnce length = " + str(config.seq_length))

    config.available_characters = [
        unambiguous_dna_alphabet,
        unambiguous_rna_alphabet,
        nucleic_alphabet,
        dna_alphabet,
        rna_alphabet,
        reduced_nucleic_alphabet,
        unambiguous_protein_alphabet,
        reduced_protein_alphabet,
        protein_alphabet,
        generic_alphabet,
    ]
    count = 0
    current_chracters = "".join(config.charaters)

    for guess in config.available_characters:
        if guess.alphabetic(current_chracters):
            config.seq_type = guess
            break
        count += 1
    if config.seq_type == "dna":
        raise Exception("No match")

    print("sequnce type = " + str(config.seq_type))
    # print(config.updated_dict['ABDA_AEDAE'].seq)

    config.seq_dict = {}  # dictionary for storing sequnce matrix
    config.seq_counter = (
        {}
    )  # dictionary storing count matrix as a whole(generate unmodifed sequnce logo)
    config.matrix = np.zeros(
        (len(config.seq_type), config.seq_length), dtype=float
    )  # a default matrix for each individual leaf/node, to store character counts.
    config.existing_characters = []

    # initialize the seq_counter dictionary
    for i in range(0, len(config.seq_type)):
        config.seq_counter[config.seq_type[i]] = np.zeros(
            (config.seq_length,), dtype=int
        ).tolist()

    start_time = time.time()
    PIC_postorder(tree)

    array = convert_count_to_array(config.seq_counter)
    print(array)

    logo_data = LogoData.from_counts(alphabet=config.seq_type, counts=array)

    logo_options = LogoOptions()
    logo_options.title = "without PIC"
    logo_options.stack_width = 50  # increase width of each position
    logo_options.stack_height = 100  # increase overall height

    logo_format = LogoFormat(logo_data, logo_options)

    # Save as PNG
    with open("Regular_logo.pdf", "wb") as g:
        g.write(pdf_formatter(logo_data, logo_format))

    # print(tree)
    # print(seq_counter)
    # counts = {'A' : [3,4,5,6], 'C': [2,3,1,1], 'T': [2,1,3,1], 'G': [3,2,1,2]}
    # print(config.seq_counter)
    """
    fig, axarr = draw_logo(config.seq_counter, data_type='counts', seq_type= config.seq_type , yaxis='probability', colorscheme='hydrophobicity')
    fig.tight_layout()
    fig.show()
    fig.savefig("without_PIC.png")
    print("without_PIC.png saved")'
    """


"""
    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(convert_count_to_scores(config.seq_counter),colorscheme='hydrophobicity',draw_axis=True)
    fig.tight_layout()
    fig.savefig("without_PIC_bits.png")
    print("without_PIC_bits.png saved")

    print("---Done in %s seconds ---" % round(time.time() - start_time))


    test = [[('A', 1.0), ('C', 0.0), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)], 
 [('A', 1.0), ('C', 0.0), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.125), ('C', 0.0), ('D', 0.125), ('E', 0.125), ('I', 0.125), ('Q', 0.0), ('R', 0.125), ('S', 0.125), ('T', 0.125), ('V', 0.125), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.125), ('C', 0.0), ('D', 0.125), ('E', 0.125), ('I', 0.125), ('Q', 0.0), ('R', 0.125), ('S', 0.125), ('T', 0.125), ('V', 0.125), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.5), ('C', 0.5), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.5), ('C', 0.5), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.125), ('C', 0.875), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.5), ('C', 0.5), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.5), ('C', 0.5), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.0), ('R', 0.0), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)],
 [('A', 0.0), ('C', 0.0), ('D', 0.0), ('E', 0.0), ('I', 0.0), ('Q', 0.125), ('R', 0.875), ('S', 0.0), ('T', 0.0), ('V', 0.0), ('N',0.0), ('G',0.0), ('H',0.0), ('L',0.0), ('K',0.0), ('M',0.0), ('F',0.0) ,('P',0.0),('W',0.0) ,('Y',0.0)]]
    
    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(test,yaxis='bits',seq_type='aa',colorscheme='chemistry',draw_axis=True)
    fig.tight_layout()
    fig.savefig("ex1_t1_new.png")
"""

if __name__ == "__main__":
    try:
        main()
    except argparse.ArgumentError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
