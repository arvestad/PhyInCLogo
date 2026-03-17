import os
import pytest
import numpy as np
from pathlib import Path
from types import SimpleNamespace

import phyinc.config as config
from phyinc.phyinc import (
    add_length,
    set_seq,
    enforce_bifurcations,
    output_info,
    validate_path,
    find_clades,
    export_scores_to_file,
    convert_matrix_to_array,
    convert_count_to_array,
)
from weblogo.seq import unambiguous_dna_alphabet

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples', 'synthetic_data')


class MockClade:
    """Minimal stand-in for a Bio.Phylo clade used in PIC calculations."""
    def __init__(self, branch_length, name=None, clades=None):
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []


# --- add_length ---

def test_add_length_both_zero():
    leaf_i = MockClade(0)
    leaf_j = MockClade(0)
    assert add_length(leaf_i, leaf_j) == 0.0


def test_add_length_equal_branches():
    leaf_i = MockClade(1.0)
    leaf_j = MockClade(1.0)
    result = add_length(leaf_i, leaf_j)
    assert result == pytest.approx(0.5)


def test_add_length_unequal_branches():
    leaf_i = MockClade(1.0)
    leaf_j = MockClade(3.0)
    # (1*3)/(1+3) = 0.75
    assert add_length(leaf_i, leaf_j) == pytest.approx(0.75)


def test_add_length_one_zero_branch():
    leaf_i = MockClade(0.0)
    leaf_j = MockClade(2.0)
    # (0*2)/(0+2) = 0
    assert add_length(leaf_i, leaf_j) == pytest.approx(0.0)


# --- set_seq ---

@pytest.fixture(autouse=True)
def reset_config():
    """Reset shared config state before each test."""
    config.seq_dict = {}
    config.matrix = None
    config.seq_type = unambiguous_dna_alphabet
    config.seq_length = 3
    yield
    config.seq_dict = {}


def test_set_seq_both_zero_returns_matrix():
    sentinel = np.zeros((4, 3))
    config.matrix = sentinel
    leaf_i = MockClade(0)
    leaf_j = MockClade(0)
    result = set_seq(leaf_i, leaf_j)
    assert result is sentinel


def test_set_seq_equal_branches():
    mat_i = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    mat_j = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    leaf_i = MockClade(1.0)
    leaf_j = MockClade(1.0)
    config.seq_dict[leaf_i] = mat_i
    config.seq_dict[leaf_j] = mat_j
    result = set_seq(leaf_i, leaf_j)
    # With equal branch lengths l=r=1: (l/(l+r))*mat_j + (r/(l+r))*mat_i = 0.5*(mat_i+mat_j)
    expected = 0.5 * (mat_i + mat_j)
    np.testing.assert_array_almost_equal(result, expected)


def test_set_seq_asymmetric_branches():
    mat_i = np.ones((4, 3))
    mat_j = np.zeros((4, 3))
    leaf_i = MockClade(3.0)
    leaf_j = MockClade(1.0)
    config.seq_dict[leaf_i] = mat_i
    config.seq_dict[leaf_j] = mat_j
    result = set_seq(leaf_i, leaf_j)
    # l=3, r=1: (l/(l+r))*mat_j + (r/(l+r))*mat_i = 0.75*zeros + 0.25*ones = 0.25
    np.testing.assert_array_almost_equal(result, 0.25 * np.ones((4, 3)))


# --- enforce_bifurcations ---

def test_enforce_bifurcations_two_children():
    mat_a = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    mat_b = np.array([[0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    config.seq_dict[child_a] = mat_a
    config.seq_dict[child_b] = mat_b
    branch, seq = enforce_bifurcations([child_a, child_b])
    assert branch == pytest.approx(0.5)
    np.testing.assert_array_almost_equal(seq, 0.5 * (mat_a + mat_b))


def test_enforce_bifurcations_both_zero_branches():
    sentinel = np.zeros((4, 2))
    config.matrix = sentinel
    child_a = MockClade(0.0)
    child_b = MockClade(0.0)
    config.seq_dict[child_a] = np.ones((4, 2))
    config.seq_dict[child_b] = np.ones((4, 2))
    branch, seq = enforce_bifurcations([child_a, child_b])
    assert branch == 0
    assert seq is sentinel


def test_enforce_bifurcations_three_children():
    mat = np.eye(4, 3)
    child_a = MockClade(1.0)
    child_b = MockClade(1.0)
    child_c = MockClade(1.0)
    config.seq_dict[child_a] = mat.copy()
    config.seq_dict[child_b] = mat.copy()
    config.seq_dict[child_c] = mat.copy()
    branch, seq = enforce_bifurcations([child_a, child_b, child_c])
    assert isinstance(branch, float)
    assert seq.shape == mat.shape


# --- convert_matrix_to_array ---

def test_convert_matrix_to_array_shape():
    config.seq_type = unambiguous_dna_alphabet  # 'ACGT', 4 chars
    config.seq_length = 5
    matrix = np.zeros((4, 5))
    result = convert_matrix_to_array(matrix)
    assert result.shape == (5, 4)


def test_convert_matrix_to_array_values():
    config.seq_type = unambiguous_dna_alphabet
    config.seq_length = 2
    matrix = np.zeros((4, 2))
    matrix[0][0] = 1.0   # A at position 0
    matrix[2][1] = 0.5   # G at position 1
    result = convert_matrix_to_array(matrix)
    assert result[0][0] == pytest.approx(1.0)   # position 0, char A (index 0)
    assert result[1][2] == pytest.approx(0.5)   # position 1, char G (index 2)


# --- convert_count_to_array ---

def test_convert_count_to_array_shape():
    config.seq_type = unambiguous_dna_alphabet
    config.seq_length = 3
    count_matrix = {
        'A': [1, 0, 2],
        'C': [0, 3, 0],
        'G': [0, 0, 1],
        'T': [2, 1, 0],
    }
    result = convert_count_to_array(count_matrix)
    assert result.shape == (3, 4)


def test_convert_count_to_array_missing_key():
    config.seq_type = unambiguous_dna_alphabet
    config.seq_length = 2
    # Only 'A' present; other characters default to 0
    count_matrix = {'A': [5, 3]}
    result = convert_count_to_array(count_matrix)
    assert result[0][0] == pytest.approx(5)   # A at pos 0
    assert result[0][1] == pytest.approx(0)   # C at pos 0 (missing -> 0)


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


# --- export_scores_to_file ---

def test_export_scores_to_file_creates_file(tmp_path):
    f = tmp_path / 'scores.txt'
    export_scores_to_file('1.23 4.56', str(f))
    assert f.read_text() == '1.23 4.56'


def test_export_scores_to_file_existing_file(tmp_path, capsys):
    f = tmp_path / 'scores.txt'
    f.write_text('old')
    export_scores_to_file('new', str(f))
    captured = capsys.readouterr()
    assert 'Already exists' in captured.out
    assert f.read_text() == 'old'  # file unchanged


# --- find_clades ---

def test_find_clades_no_match():
    root = MockClade(0, name='root', clades=[
        MockClade(1.0, name='A'),
        MockClade(1.0, name='B'),
    ])
    result = find_clades(root, lambda c: c.name == 'Z')
    assert result == []


def test_find_clades_matches_root():
    root = MockClade(0, name='root')
    result = find_clades(root, lambda c: c.name == 'root')
    assert len(result) == 1
    assert result[0] is root


def test_find_clades_matches_leaf():
    leaf = MockClade(1.0, name='leaf')
    root = MockClade(0, name='root', clades=[leaf])
    result = find_clades(root, lambda c: c.name == 'leaf')
    assert len(result) == 1
    assert result[0] is leaf


def test_find_clades_multiple_matches():
    leaf_a = MockClade(1.0, name='leaf')
    leaf_b = MockClade(2.0, name='leaf')
    root = MockClade(0, name='root', clades=[leaf_a, leaf_b])
    result = find_clades(root, lambda c: c.name == 'leaf')
    assert len(result) == 2


def test_find_clades_nested():
    grandchild = MockClade(0.5, name='target')
    child = MockClade(1.0, name='child', clades=[grandchild])
    root = MockClade(0, name='root', clades=[child])
    result = find_clades(root, lambda c: c.name == 'target')
    assert len(result) == 1
    assert result[0] is grandchild
