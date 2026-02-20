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



def read_sequences(filename, filetype="fasta", args) -> tuple[dict, Alphabet]:
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


def infer_sequence_type(char_set, args):

    seq_type = 'dna'           # The default
    # 1. Guess type by comparing the character set with WebLogo alphabets

    # 2. If user desires another sequence type, then ensure seqs are compatible
    if args.type:
        pass
    else:
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
    current_chracters = "".join(config.characters)

    for guess in config.available_characters:
        if guess.alphabetic(current_chracters):
            config.seq_type = guess
            break
    if config.seq_type == "dna":
        raise Exception("No match")
    
