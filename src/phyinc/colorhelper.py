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

# The Wes Anderson palette is borrowed from the R package 'wesanderson'
# with the help from Claude Code.
wes_anderson = {                
    'A': "#02401B",  # 16 Jungle bottle      — Darjeeling Ltd.
    'C': "#7294D4",  # 5  Concierge blue     — Grand Budapest
    'D': "#FD6467",  # 1  Futura red         — Grand Budapest
    'E': "#F1BB7B",  # 2  Mendl's peach      — Grand Budapest
    'F': "#5B1A18",  # 3  Lobby maroon       — Grand Budapest
    'G': "#D8B70A",  # 15 Saffron road       — Darjeeling Ltd.
    'H': "#3B9AB2",  # 6  Zissou sea         — Life Aquatic
    'I': "#35274A",  # 7  Blume purple       — Rushmore
    'K': "#F21A00",  # 8  Jaguar orange      — Life Aquatic
    'L': "#F3DF6C",  # 9  Moonrise gold      — Moonrise Kingdom
    'M': "#85D4E3",  # 10 Island sky         — Moonrise Kingdom
    'N': "#9C964A",  # 11 Tent khaki         — Moonrise Kingdom
    'P': "#899DA4",  # 12 Tenenbaum slate    — Royal Tenenbaums
    'Q': "#A2A475",  # 17 Dusty sage         — Darjeeling Ltd.
    'R': "#DC863B",  # 14 Etheline amber     — Royal Tenenbaums
    'S': "#E6A0C4",  # 4  Agatha pink        — Grand Budapest
    'T': "#C93312",  # 13 Richie's headband  — Royal Tenenbaums
    'V': "#46ACC8",  # 18 Fox sky cyan       — Fantastic Mr. Fox
    'W': "#6C6F39",  # 19 Autumn olive       — Fantastic Mr. Fox
    'Y': "#446455",  # 20 Chevalier pine     — The Chevalier
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
    'wesanderson':    wes_anderson,
    'random':         'random',
}


def decide_color_scheme(args, color_scheme):
    if args.color_scheme == 'guess':
        return color_scheme
    else:
        return color_lookup.get(args.color_scheme, 'gray')
