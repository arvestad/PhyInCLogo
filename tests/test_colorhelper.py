import pytest
from types import SimpleNamespace

from phyinc.colorhelper import decide_color_scheme
from weblogo.colorscheme import (
    monochrome,
    nucleotide,
    hydrophobicity,
    chemistry,
    charge,
    taylor,
)


def test_decide_color_scheme_guess_returns_passed_scheme():
    args = SimpleNamespace(color_scheme='guess')
    result = decide_color_scheme(args, nucleotide)
    assert result is nucleotide


def test_decide_color_scheme_guess_preserves_hydrophobicity():
    args = SimpleNamespace(color_scheme='guess')
    result = decide_color_scheme(args, hydrophobicity)
    assert result is hydrophobicity


def test_decide_color_scheme_monochrome():
    args = SimpleNamespace(color_scheme='monochrome')
    result = decide_color_scheme(args, nucleotide)
    assert result is monochrome


def test_decide_color_scheme_nucleotide():
    args = SimpleNamespace(color_scheme='nucleotide')
    result = decide_color_scheme(args, hydrophobicity)
    assert result is nucleotide


def test_decide_color_scheme_hydrophobicity():
    args = SimpleNamespace(color_scheme='hydrophobicity')
    result = decide_color_scheme(args, nucleotide)
    assert result is hydrophobicity


def test_decide_color_scheme_chemistry():
    args = SimpleNamespace(color_scheme='chemistry')
    result = decide_color_scheme(args, nucleotide)
    assert result is chemistry


def test_decide_color_scheme_charge():
    args = SimpleNamespace(color_scheme='charge')
    result = decide_color_scheme(args, nucleotide)
    assert result is charge


def test_decide_color_scheme_taylor():
    args = SimpleNamespace(color_scheme='taylor')
    result = decide_color_scheme(args, nucleotide)
    assert result is taylor


def test_decide_color_scheme_unknown_falls_back_to_monochrome():
    args = SimpleNamespace(color_scheme='nonexistent')
    result = decide_color_scheme(args, nucleotide)
    assert result == monochrome
