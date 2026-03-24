import argparse
import logging
import sys

import numpy as np

import matplotlib
matplotlib.use("TkAgg")

from pathlib import Path
from Bio import Phylo

import phyinc.io as io
from phyinc.colorhelper import decide_color_scheme
from importlib.metadata import version

from weblogo import LogoData, LogoOptions, LogoFormat, formatters as logo_formatters


output_formats = ["pdf", "eps", "png", "jpeg", "gif"]
fineprint = "Created with PhyInC Logo + WebLogo"
description_str = """
Make sequence logos using Felsenstain's phylogenetically independent contrast metod to take evolution into account.
"""

def setup_argparse():
    # Apply PIC on tree file and corresbond Name-sequnce file, outputs 2 .png file(sequnce logo with and with out applying PIC).
    parser = argparse.ArgumentParser(
        description=description_str + f" Version {version('phyinc')}."
    )
    parser.add_argument(
        "seq_filename",
        type=str,
        help="Path to the Fasta file. Assumes sequence name in the format of: ETA_STAAU/96-110",
    )
    parser.add_argument(
        "tree_filename", type=str, help="Path to the .tree file (newick format)"
    )
    parser.add_argument(
        "-c",
        "--color-scheme",
        choices=["monochrome", "nucleotide", "base_pairing", "hydrophobicity", "chemistry", "charge", "taylor", "guess"],
        default="guess",
        help="Choose color scheme. If using 'guess', then sequence type chooses between 'nucleotide' and 'taylor'."
    )
    parser.add_argument(
        "--coords",
        action="store_true",
        help="Ignore domain coordinates in the sequence file when matching accessions to the treefile. I.e., for an accession like 'ACC/17-33' only 'ACC' is used."
    )
    parser.add_argument(
        "-t",
        "--type",
        choices=["aa", "dna", "rna", "guess"],
        default="guess",
        help="Specify what sequence type to assume. Default: %(default)s")
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Name of outfile. If this option is not used, it will be inferred from the sequence file. If the filename ends with .pdf, .png, or corresponding to another accepted format, that output format will be chosen.")
    parser.add_argument(
        "-f",
        "--format",
        choices=output_formats,
        default="pdf",
        help=f"Choose an output format, one of {output_formats}. This option is ignored if --outfile is used and a format is given in the filename.",
    )
    parser.add_argument(
        "-n",
        "--no-fineprint",
        action="store_true",
        help=f'Do not add a string indicating what software produced the logo ("{fineprint}")',
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=f"Write more information to stderr",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"phyinc version {version('phyinc')}"
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


def convert_matrix_to_array(matrix, seq_type, seq_length):
    array = []
    for i in range(0, seq_length):
        sub_array = []
        for c in seq_type:
            sub_array.append(matrix[seq_type.letters().index(c)][i])
        array.append(sub_array)
    return np.array(array)



def set_seq(leaf_i, leaf_j, seq_dict, zero_matrix):
    """Calculate frequency matrix (x_i) for parent node"""
    # TODO: choose better function name
    l = float(leaf_i.branch_length)
    r = float(leaf_j.branch_length)
    if (l == 0) and (r == 0):
        return zero_matrix
    else:
        seq_matrix = np.add(
            (l / (l + r)) * seq_dict[leaf_j],
            (r / (l + r)) * seq_dict[leaf_i],
        )
        return seq_matrix


def add_length(leaf_i, leaf_j):
    """calculate new branch length (v_i') for parent node"""
    l = float(leaf_i.branch_length)
    r = float(leaf_j.branch_length)
    if (l == 0) and (r == 0):
        return float(0)
    else:
        return (l * r) / (l + r)


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


def pic_seqlogo(tree, logo_formatter, color_scheme, args, alignment, seq_type, seq_length):
    # Build local state for this run — no global config needed.
    zero_matrix = np.zeros((len(seq_type), seq_length), dtype=float)
    seq_dict = {}

    for child in tree.clade:
        traverse_postorder(child, alignment, seq_type, zero_matrix, seq_dict)

    result, seq_matrix = collapse_bifurcations(tree.clade, seq_dict, zero_matrix)
    array = convert_matrix_to_array(seq_matrix, seq_type, seq_length)

    logo_data = LogoData.from_counts(alphabet=seq_type, counts=array)
    logo_options = LogoOptions()

    if args.no_fineprint:
        logo_options.show_fineprint = False
    else:
        # add the fine print
        logo_options.fineprint = fineprint

    logo_options.color_scheme = color_scheme
    logo_options.stack_width = 50  # increase width of each position
    logo_options.stack_height = 100  # increase overall height

    logo_format = LogoFormat(logo_data, logo_options)
    return logo_formatter(logo_data, logo_format)


def traverse_postorder(clade, alignment, seq_type, zero_matrix, seq_dict):
    if len(clade) == 0:  # only tips of the tree will have length 0
        clade.seq = str(
            alignment[clade.name].seq
        )  # store str(sequences) as an attribute for the clade object (of biopython package).

        seq_matrix = zero_matrix.copy()  # make a copy of the default sequence matrix
        for i in range(0, len(clade.seq)):
            character_index = seq_type.letters().index(clade.seq[i].upper())
            seq_matrix[character_index, i] = float(
                1
            )  # a count matrix for each individual leaf, used later to calculate frequency matrix for parent nodes

        seq_dict[clade] = seq_matrix  # stores the matrix to dictionary
    else:  # a parent node, because len(clade) > 0
        for child in clade:
            traverse_postorder(child, alignment, seq_type, zero_matrix, seq_dict)
        if len(clade) == 2:
            clade.branch_length = float(clade.branch_length) + add_length(
                clade[0], clade[1]
            )
            seq_dict[clade] = set_seq(clade[0], clade[1], seq_dict, zero_matrix)
        else:
            branch_length, seq_matrix = collapse_bifurcations(clade, seq_dict, zero_matrix)
            clade.branch_length = float(clade.branch_length) + branch_length
            seq_dict[clade] = seq_matrix


def validate_path(filename):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"Cannot open {filename}.")
    return path


def main():
    # Set up argument parser
    ap = setup_argparse()
    args = ap.parse_args()

    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format="phyinc: %(message)s",
        )

    try:
        tree_file = validate_path(args.tree_filename)
        seq_file = validate_path(args.seq_filename)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    outfilename, logo_artist = output_info(args)  # Infer details before computing

    try:
        tree = Phylo.read(tree_file, "newick")
    except IOError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(5)

    try:
        alignment, seq_type, seq_length, color_scheme = io.read_sequences(seq_file, "fasta", args)
    except IOError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(3)
    except KeyError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(4)

    logging.info(f"Alignment width:  {seq_length}")
    logging.info(f"Sequence type:    {seq_type}")
    logging.info(f"Number of leaves: {tree.count_terminals()}")

    try:
        io.check_accession_consistency(alignment, tree, args.coords)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(4)

    color_scheme = decide_color_scheme(args, color_scheme)

    try:
        logo = pic_seqlogo(tree, logo_artist, color_scheme, args,
                           alignment, seq_type, seq_length)
        with open(outfilename, "wb") as out:
            out.write(logo)
    except KeyError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
        


if __name__ == "__main__":
    main()
        
