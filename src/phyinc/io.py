from Bio import SeqIO
import phyinc.config as config

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
nucleotide_alphabets = [
        unambiguous_dna_alphabet,
        unambiguous_rna_alphabet,
        reduced_nucleic_alphabet,
        dna_alphabet,
        rna_alphabet,
        nucleic_alphabet,
]
protein_alphabets = [
        unambiguous_protein_alphabet,
        reduced_protein_alphabet,
        protein_alphabet,
]

def infer_sequence_type(char_set, args):

    # 1. Guess type by comparing the character set with WebLogo alphabets
    seq_type = None
    available_alphabets = nucleotide_alphabets + protein_alphabets
    for alphabet in available_alphabets:
        if char_set.issubset(set(alphabet)):
            seq_type = alphabet
            break               # Found it!

    # 2. If user desires another sequence type, then ensure seqs are compatible
    if args.type:
        if args.type == 'dna':
            pass
        elif args.type == 'rna':
            pass
        elif args.type == 'aa':
            pass
        else:
            pass                # We end up here if == 'guess'

    return seq_type


    # Old version:
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
    for sigma in config.available_characters:
        print(f'{type(sigma)}\t {sigma}')

    current_characters = "".join(config.characters)

    for guess in config.available_characters:
        if guess.alphabetic(current_characters):
            config.seq_type = guess
            break
    if config.seq_type == "dna":
        raise Exception("No match")
    

def read_sequences(filename, filetype, args):
    """
    Read the input alignment and perform basic length checks.
    Also, determine what kind of bio sequence we read.
    Returns a dict with the alignment and """

    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    seq_dict = {}
    characters = set()

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

        characters.update(set(value.upper()))
        seq_dict[key] = value

    seq_type = infer_sequence_type(characters, args)
    config.seq_type = seq_type     # TODO: remove!
    config.characters = characters # TODO: remove!
    return seq_dict, seq_type


