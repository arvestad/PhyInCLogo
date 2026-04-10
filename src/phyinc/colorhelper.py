import random
import matplotlib.cm as cm


def random_color_scheme(letters):
    """
    Map each character in letters to a distinct color drawn from tab20.
    The assignment is shuffled so repeated calls give different palettes.
    """
    cmap = cm.get_cmap('tab20', 20)
    colors = [cmap(i) for i in range(20)]
    random.shuffle(colors)
    return {
        char: '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))
        for char, (r, g, b, _) in zip(letters, colors)
    }


# Taylor (1997) amino acid color scheme
taylor = {
    'A': '#CCFF00', 'C': '#FFFF00', 'D': '#FF0000', 'E': '#FF0066',
    'F': '#00FF66', 'G': '#FF9900', 'H': '#0066FF', 'I': '#66FF00',
    'K': '#6600FF', 'L': '#33FF00', 'M': '#00FF00', 'N': '#CC00FF',
    'P': '#FFAA00', 'Q': '#FF00CC', 'R': '#0000FF', 'S': '#FF3300',
    'T': '#FF6600', 'V': '#99FF00', 'W': '#00CCFF', 'Y': '#00FFCC',
}

# Maps --color-scheme argument values to logomaker color_scheme values.
# Strings are logomaker built-ins; dicts are character-to-color mappings.
color_lookup = {
    'gray':           'gray',
    'nucleotide':     'classic',
    'base_pairing':   'base_pairing',
    'hydrophobicity': 'hydrophobicity',
    'chemistry':      'chemistry',
    'charge':         'charge',
    'taylor':         taylor,
    'random':         'random',
}


def decide_color_scheme(args, color_scheme):
    if args.color_scheme == 'guess':
        return color_scheme
    else:
        return color_lookup.get(args.color_scheme, 'gray')
