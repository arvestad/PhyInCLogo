from Bio import SeqIO
import logging
import phyinc.config as config

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
        choice = match_alphabets(char_set, unambiguous_dna_alphabet)
        if choice is not None:
            return choice, 'nucleotide'
        choice = match_alphabets(char_set, unambiguous_rna_alphabet)
        if choice is not None:
            return choice, 'nucleotide'
        return match_alphabets(char_set, protein_alphabets), hydrophobicity

    elif args.type == 'dna':
        return match_alphabets(char_set, dna_alphabets), nucleotide
    elif args.type == 'rna':
        return match_alphabets(char_set, rna_alphabets), nucleotide
    elif args.type == 'aa':
        return match_alphabets(char_set, aa_alphabets), hydrophobicity
    else:
        raise Exception("Bug in infer_sequence_type. Please report!")
    

def read_sequences(filename, filetype, args):
    """
    Read the input alignment and perform basic length checks.
    Also, determine what kind of bio sequence we read.

    Returns a dict with the alignment, a matching alphabet, and
    a suggested coloring scheme.
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    seq_dict = {}
    characters = set()

    alignment_width = None
    for acc, seq_str in record_dict.items():
        if len(seq_str) != alignment_width:
            if alignment_width:
                raise Exception(
                    "Sequence lengths from provided fastafile are inconsistent"
                )
            else:
                alignment_width = len(seq_str)
                config.seq_length = len(seq_str)

        seq_str = seq_str.upper()
        characters.update(set(seq_str))
        seq_dict[acc] = seq_str

    seq_type, coloring_scheme = infer_sequence_type(characters, args)
    config.seq_type = seq_type     # TODO: remove!
    config.characters = characters # TODO: remove!
    return seq_dict, seq_type, coloring_scheme


