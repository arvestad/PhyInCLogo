from phyinc.io import infer_sequence_type, match_alphabets

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


def test_dna_inference():
    alpha = match_alphabets(set('ACGT'), dna_alphabets)
    assert str(alpha) == 'ACGT'

    for c in 'ACGT': 
        alpha = match_alphabets(set(c), dna_alphabets)
        assert str(alpha) == 'ACGT'
        
