from weblogo.colorscheme import (
    monochrome,
    nucleotide,
    base_pairing,
    hydrophobicity,
    chemistry,
    charge,
    taylor
    )

color_lookup = {
    'monochrome': monochrome,
    'nucleotide': nucleotide,
    'base_pairing': base_pairing,
    'hydrophobicity': hydrophobicity,
    'chemistry': chemistry,
    'charge': charge,
    'taylor': taylor,
    }

def decide_color_scheme(args, color_scheme):
    if args.color_scheme == 'guess':
        return color_scheme
    else:
        return color_lookup.get(args.color_scheme, 'monochrome')

