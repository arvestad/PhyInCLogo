import argparse
import config
import logging
import math
import os
import sys
import weblogo  # is this needed?

import numpy as np
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from pathlib import Path
from Bio import Phylo
from Bio import SeqIO


from weblogo import LogoData, LogoOptions, LogoFormat, formatters as logo_formatters
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


output_formats = ["pdf", "eps", "png", "jpeg", "gif"]
fineprint = "Created with PhyInC Logo + WebLogo"


def setup_argparse():
    # Apply PIC on tree file and corresbond Name-sequnce file, outputs 2 .png file(sequnce logo with and with out applying PIC).
    parser = argparse.ArgumentParser(
        description="Make sequence logos using Felsenstain's phylogenetically independent contrast metod to take evolution into account."
    )
    parser.add_argument(
        "tree_filename", type=str, help="Path to the .tree file (newick format)"
    )
    parser.add_argument(
        "seq_filename",
        type=str,
        help="Path to the Fasta file. Assumes sequence name in the format of: ETA_STAAU/96-110",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Name of outfile. If this option is not used, it will be inferred from the sequence file. If the filename ends with .pdf, .png, or corresponding to another accepted format, that output format will be chosen.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=output_formats,
        help=f"Choose an output format, one of {output_formats}. This option is ignored if --outfile is used and a format is given in the filename.",
    )
    parser.add_argument(
        "-n",
        "--no-fineprint",
        action="store_true",
        help=f'Do not add a string indicating what software produced the logo ("{fineprint}")',
    )

    return parser


def output_info(args):
    """
    Figure out what file to write logo to and what format to use.
    Returns the pair `filename`, `formatter` and ensures they make sense together.
    The formatter converts logo data to a picture.
    """
    if args.outfile:
        outfile = Path(args.outfile)
        graphics_format = outfile.suffix.strip(".")
        if not graphics_format:
            raise ValueError(
                f"Give your outfile a suffix that determines output format, one of {output_formats}"
            )
        if not graphics_format in output_formats:
            raise ValueError(
                f'"{graphics_format}" is not a valid output format. Use one of {output_formats}.'
            )
    else:
        graphics_format = "pdf"  # Good default
        if args.format:
            graphics_format = (
                args.format
            )  # If argparse accepted the string, then it is good

        outfile = args.seq_filename + "_logo." + graphics_format
    return outfile, logo_formatters[graphics_format]


def export_scores_to_file(scores_str, file_name):
    try:
        with open(file_name, "x") as f:
            f.write(scores_str)
    except FileExistsError:
        print("Already exists.")


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


def enforce_bifurcations(children):
    """
    In case the tree has multifurcating nodes, we create bifurcations
    by adding very short edges.
    """
    branch_length_temp = []
    seq_matrix_temp = []
    for child in children:
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


def pic_seqlogo(tree, logo_formatter, no_fine_print):

    for child in tree.clade:
        traverse_postorder(child)

    result, seq_matrix = enforce_bifurcations(tree.clade)

    array = convert_matrix_to_array(seq_matrix)

    logo_data = LogoData.from_counts(alphabet=config.seq_type, counts=array)

    logo_options = LogoOptions()
    # logo_options.title = "With PIC logo"

    if no_fine_print:
        logo_options.show_fineprint = False

    else:
        # add the fine print
        logo_options.fineprint = fineprint

    logo_options.stack_width = 50  # increase width of each position
    logo_options.stack_height = 100  # increase overall height

    logo_format = LogoFormat(logo_data, logo_options)
    return logo_formatter(logo_data, logo_format)

    # Save as PNG
    # with open("With_PIC_logo.png", "wb") as f:
    #     f.write()


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
            branch_length, seq_matrix = enforce_bifurcations(clade)
            clade.branch_length = float(clade.branch_length) + branch_length
            config.seq_dict[clade] = seq_matrix


def read_sequences(filename, filetype="fasta"):
    """checks length and type of sequnces."""

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    seq_dict = {}

    alignment_width = None
    for key, value in record_dict.items():
        if len(value) != alignment_width:
            if alignment_width:
                raise Exception(
                    "Sequence length from provided fastafile are inconsistent"
                )
            else:
                alignment_width = len(value)
                config.seq_length = len(value)

        config.characters.update(set(value.upper()))
        seq_dict[key] = value

    return seq_dict


def find_clades(clade, condition):
    """Find clades matching a condition, used in testing."""
    matches = []
    if condition(clade):
        matches.append(clade)
    for sub_clade in clade.clades:
        matches.extend(find_clades(sub_clade, condition))
    return matches


def validate_path(filename):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f".tree {filename} is not a file.")
    return path


def main():
    # Set up argument parser
    ap = setup_argparse()
    args = ap.parse_args()

    tree_file = validate_path(args.tree_filename)
    seq_file = validate_path(args.seq_filename)

    outfilename, logo_artist = output_info(args)  # Infer details before computing

    tree = Phylo.read(tree_file, "newick")
    config.terminals = tree.count_terminals()

    # config.py to store global variables
    config.updated_dict = read_sequences(seq_file, "fasta")
    logging.info("Alignment width = " + str(config.seq_length))

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
    current_chracters = "".join(config.characters)

    for guess in config.available_characters:
        if guess.alphabetic(current_chracters):
            config.seq_type = guess
            break
        count += 1
    if config.seq_type == "dna":
        raise Exception("No match")

    logging.info("Sequence type = " + str(config.seq_type))

    config.seq_dict = {}  # dictionary for storing sequnce matrix
    config.seq_counter = (
        {}
    )  # dictionary storing count matrix as a whole(generate unmodifed sequnce logo)
    config.matrix = np.zeros(
        (len(config.seq_type), config.seq_length), dtype=float
    )  # a default matrix for each individual leaf/node, to store character counts.
    config.existing_characters = []

    # initialize the seq_counter dictionary
    # length_of_config = len(config.seq_type.letters())
    # print(length_of_config)
    # print(type(config.seq_type), config.seq_type)
    # print(config.seq_type.letters()[1])

    for i in range(0, len(config.seq_type.letters())):
        # Alphabet object is not subscriptable so you cannot call config.seq_type[i]
        config.seq_counter[config.seq_type.letters()[i]] = np.zeros(
            (config.seq_length,), dtype=int
        ).tolist()

    logo = pic_seqlogo(tree, logo_artist, args.no_fineprint)
    with open(outfilename, "wb") as out:
        out.write(logo)

    # array = convert_count_to_array(config.seq_counter)

    # logo_data = LogoData.from_counts(alphabet=config.seq_type, counts=array)

    # logo_options = LogoOptions()
    # logo_options.title = "without PIC"
    # logo_options.stack_width = 50  # increase width of each position
    # logo_options.stack_height = 100  # increase overall height

    # logo_format = LogoFormat(logo_data, logo_options)

    # Save as PNG
    # with open("Regular_logo.pdf", "wb") as g:
    #     g.write(pdf_formatter(logo_data, logo_format))


if __name__ == "__main__":
    try:
        main()
    except argparse.ArgumentError as e:
        print(f"Error: {e}", file=sys.stderr)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
