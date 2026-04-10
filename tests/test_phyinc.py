import os
import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace

from phyinc.phyinc import (
    collapse_bifurcations,
    compute_pic_array,
    convert_matrix_to_array,
)
from phyinc.main import output_info, validate_path
from phyinc.io import unambiguous_dna_alphabet

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples', 'synthetic_data')


class MockClade:
    """Minimal stand-in for a Bio.Phylo clade used in PIC calculations."""
    def __init__(self, branch_length, name=None, clades=None):
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []



# --- collapse_bifurcations ---

def test_collapse_bifurcations_two_children():
    mat_a = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    mat_b = np.array([[0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    seq_dict = {child_a: mat_a, child_b: mat_b}
    branch, seq = collapse_bifurcations([child_a, child_b], seq_dict, np.zeros((4, 2)))
    assert branch == pytest.approx(0.5)
    np.testing.assert_array_almost_equal(seq, 0.5 * (mat_a + mat_b))


def test_collapse_bifurcations_both_zero_branches():
    sentinel = np.zeros((4, 2))
    child_a = MockClade(0.0)
    child_b = MockClade(0.0)
    seq_dict = {child_a: np.ones((4, 2)), child_b: np.ones((4, 2))}
    branch, seq = collapse_bifurcations([child_a, child_b], seq_dict, sentinel)
    assert branch == 0
    assert seq is sentinel


def test_collapse_bifurcations_three_children():
    mat = np.eye(4, 3)
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    child_c = MockClade(1.0)
    seq_dict = {child_a: mat.copy(), child_b: mat.copy(), child_c: mat.copy()}
    branch, seq = collapse_bifurcations([child_a, child_b, child_c], seq_dict, np.zeros((4, 3)))
    assert isinstance(branch, float)
    assert seq.shape == mat.shape


# --- convert_matrix_to_array ---

def test_convert_matrix_to_array_shape():
    matrix = np.zeros((4, 5))
    result = convert_matrix_to_array(matrix, unambiguous_dna_alphabet, 5)
    assert result.shape == (5, 4)


def test_convert_matrix_to_array_values():
    matrix = np.zeros((4, 2))
    matrix[0][0] = 1.0   # A at position 0
    matrix[2][1] = 0.5   # G at position 1
    result = convert_matrix_to_array(matrix, unambiguous_dna_alphabet, 2)
    assert result[0][0] == pytest.approx(1.0)   # position 0, char A (index 0)
    assert result[1][2] == pytest.approx(0.5)   # position 1, char G (index 2)



# --- output_info ---

def test_output_info_default_format():
    args = SimpleNamespace(outfile=None, format=None, seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert str(outfile) == 'seqs.fa_logo.pdf'
    assert formatter is not None


def test_output_info_explicit_format():
    args = SimpleNamespace(outfile=None, format='png', seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert str(outfile) == 'seqs.fa_logo.png'


def test_output_info_explicit_outfile():
    args = SimpleNamespace(outfile='result.eps', format=None, seq_filename='seqs.fa')
    outfile, formatter = output_info(args)
    assert Path(outfile).name == 'result.eps'


def test_output_info_outfile_no_suffix():
    args = SimpleNamespace(outfile='result', format=None, seq_filename='seqs.fa')
    with pytest.raises(ValueError, match="suffix"):
        output_info(args)


def test_output_info_outfile_unsupported_format():
    args = SimpleNamespace(outfile='result.bmp', format=None, seq_filename='seqs.fa')
    with pytest.raises(ValueError, match="not a valid output format"):
        output_info(args)


# --- validate_path ---

def test_validate_path_existing_file(tmp_path):
    f = tmp_path / 'test.txt'
    f.write_text('hello')
    result = validate_path(str(f))
    assert result == Path(str(f))


def test_validate_path_nonexistent():
    with pytest.raises(FileNotFoundError):
        validate_path('/no/such/file.tree')


def test_validate_path_directory(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_path(str(tmp_path))


# --- compute_pic_array ---

def test_compute_pic_array_shape():
    from types import SimpleNamespace
    from Bio import Phylo, SeqIO
    import io as _io
    import phyinc.io as phyinc_io

    fa_file = os.path.join(EXAMPLES_DIR, 'ex1.fa')
    tree_file = os.path.join(EXAMPLES_DIR, 'ex1_t1.tree')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    alignment, seq_type, seq_length, _ = phyinc_io.read_sequences(fa_file, 'fasta', args)
    tree = Phylo.read(tree_file, 'newick')

    array = compute_pic_array(tree, alignment, seq_type, seq_length)

    assert array.shape == (seq_length, len(seq_type))


def test_compute_pic_array_values_sum_to_one():
    from types import SimpleNamespace
    from Bio import Phylo
    import phyinc.io as phyinc_io

    fa_file = os.path.join(EXAMPLES_DIR, 'ex1.fa')
    tree_file = os.path.join(EXAMPLES_DIR, 'ex1_t1.tree')
    args = SimpleNamespace(type='aa', ignore_coords=False)
    alignment, seq_type, seq_length, _ = phyinc_io.read_sequences(fa_file, 'fasta', args)
    tree = Phylo.read(tree_file, 'newick')

    array = compute_pic_array(tree, alignment, seq_type, seq_length)

    # Each row is a weighted average of one-hot vectors, so rows should sum to 1.
    row_sums = array.sum(axis=1)
    np.testing.assert_allclose(row_sums, np.ones(seq_length), atol=1e-6)




