from Bio import SeqIO
import logging
import re

from weblogo.colorscheme import hydrophobicity, nucleotide # Default color schemes

from weblogo.seq import (
#    generic_alphabet, # I do not want to use the generic alphabet — it allows _everything_
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
dna_alphabets = [
        unambiguous_dna_alphabet,
        reduced_nucleic_alphabet,
        dna_alphabet,
        nucleic_alphabet,
]
rna_alphabets = [
        unambiguous_rna_alphabet,
        rna_alphabet,
]
protein_alphabets = [
        unambiguous_protein_alphabet,
        reduced_protein_alphabet,
        protein_alphabet,
]

all_alphabets = dna_alphabets + rna_alphabets + protein_alphabets

domain_coord_pattern = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)/(\d+)-(\d+)$')


def match_alphabets(char_set, available_alphabets):
    for alphabet in available_alphabets:
        if char_set.issubset(set(alphabet)):
            return alphabet     # Found it!
    return None


def infer_sequence_type(char_set, args):
    """
    Guess sequence type by matching alphabets.
    If the user gave a command-line instruction, then
    choose an alphabet that best matches the observed
    character set. If the input contains ambiguity
    symbols, we will assume it is a protein sequence
    and the user will need to specify the sequence type.

    We will also suggest a coloring scheme.

    This is most likely not an ideal solution but will do for now.

    Returns: an weblogo.Alphabet object and a coloring scheme.
    """
    if args.type == 'guess':           # User leaves it to us to choose
        choice = match_alphabets(char_set, all_alphabets)
        if choice is not None:
            return choice, nucleotide
        choice = match_alphabets(char_set, unambiguous_rna_alphabet)
        if choice is not None:
            return choice, nucleotide
        return match_alphabets(char_set, protein_alphabets), hydrophobicity

    elif args.type == 'dna':
        return match_alphabets(char_set, dna_alphabets), nucleotide
    elif args.type == 'rna':
        return match_alphabets(char_set, rna_alphabets), nucleotide
    elif args.type == 'aa':
        return match_alphabets(char_set, protein_alphabets), hydrophobicity
    else:
        raise Exception("Bug in infer_sequence_type. Please report!")
    

def read_sequences(filename, filetype, args):
    """
    Read the input alignment and perform basic length checks.
    Also, determine what kind of bio sequence we read.

    Returns seq_dict, seq_type, seq_length, and a suggested coloring scheme.
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    seq_dict = {}
    characters = set()

    seq_length = None
    for acc, seq_str in record_dict.items():
        if len(seq_str) != seq_length:
            if seq_length:
                raise Exception(
                    "Sequence lengths from provided fastafile are inconsistent"
                )
            else:
                seq_length = len(seq_str)

        seq_str = seq_str.upper()
        characters.update(set(seq_str))
        if args.coords:
            m = domain_coord_pattern.match(acc) # Is acc on the form "HUBBA/17-35"?
            if m:
                prot_acc, domain_start, domain_end = m.groups()
                if prot_acc in seq_dict:
                    raise IOError(f"'{prot_acc}' is a protein appearing twice, probably because you have two domains from the same protein in the input. If so, you must submit a tree inferred on the domain sequences, not on the proteins.")
                seq_dict[prot_acc] = seq_str # Map protein accession to seq

        # Always insert for full accession, even if domain start and end are included
        seq_dict[acc] = seq_str      # Map domain accession to seq

    seq_type, coloring_scheme = infer_sequence_type(characters, args)
    return seq_dict, seq_type, seq_length, coloring_scheme


def check_accession_consistency(seq_dict, tree, ignore_domain_coords=False):
    """
    Verify that sequence accessions and tree leaf names agree.

    Raises ValueError if any sequence accession has no corresponding
    leaf in the tree.
    """
    tree_leaves = {clade.name for clade in tree.get_terminals()}
    seq_accessions = set(seq_dict.keys())
    if ignore_domain_coords:
        for acc in seq_dict.keys():
            m = domain_coord_pattern.match(acc)
            if m:
                seq_accessions.remove(acc)

    missing_from_seq = tree_leaves - seq_accessions
    missing_from_tree = seq_accessions - tree_leaves

    messages = []
    if missing_from_seq:
        n = len(missing_from_seq)
        if n > 4:
            name_lst = sorted(list(missing_from_seq))
            names = ", ".join(name_lst[:4])
            messages.append(f"{n} tree leaves not found in sequence file: {names} and {n-3} more.")
        else:
            names = ", ".join(sorted(missing_from_seq))
            messages.append(f"Tree leaves not found in sequence file: {names}")

    if missing_from_tree:
        n = len(missing_from_tree)
        if n > 4:
            name_lst = sorted(list(missing_from_tree))
            names = ", ".join(name_lst[:4])
            messages.append(f"{n} sequence accessions not found in tree: {names} and {n-3} more.")
        else:
            name_lst = sorted(list(missing_from_tree))
            names = ", ".join(name_lst[:4])
            messages.append(f"Sequence accessions not found in tree: {names}")

    if messages:
        raise ValueError("\n".join(messages))

