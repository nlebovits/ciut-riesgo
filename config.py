"""
Visualization configuration for CIUT Riesgo project.
Contains color palettes, basemap settings, and line style definitions.
"""

import contextily as ctx

# Continuous color palettes
# Blue palette (5 sequential shades)
BLUE_PALETTE = [
    '#0b2c75',  # Darkest blue
    '#4d5496',  # Dark blue
    '#7f7fb8',  # Medium-dark blue
    '#afaedb',  # Medium blue
    '#e1dfff'   # Light blue
]

# Pink palette (5 sequential shades)
PINK_PALETTE = [
    '#750b40',  # Darkest pink
    '#984466',  # Dark pink
    '#bb748f',  # Medium-dark pink
    '#dda4b9',  # Medium pink
    '#ffd5e5'   # Light pink
]

# Divergent palette: pink -> neutral -> blue
DIVERGENT_PALETTE = [
    '#750b40',  # Darkest pink
    '#a25e76',  # Dark pink
    '#cba6b1',  # Light pink
    '#f1f1f1',  # Neutral
    '#aaa9c7',  # Light blue
    '#65679e',  # Medium blue
    '#0b2c75'   # Darkest blue
]

# Categorical color scheme (6 colors)
# Blue and pink from above + 4 colors that complement the scheme
CATEGORICAL_COLORS = [
    '#0b2c75',  # Blue alta
    '#750b40',  # Pink alta
    '#2d8a5c',  # Teal (complements both blue and pink)
    '#d17d2d',  # Warm orange (contrasts well)
    '#2d8a8a',  # Cyan (cool tone, complements the palette)
    '#5c2d8a'   # Purple (bridges blue and pink)
]

# Basemap configuration
BASEMAP = ctx.providers.CartoDB.PositronNoLabels

# Line weights and styles (to be fleshed out)
LINE_WEIGHTS = {
    'thin': 0.5,
    'normal': 1.0,
    'thick': 2.0,
    'very_thick': 3.0
}

LINE_STYLES = {
    'solid': '-',
    'dashed': '--',
    'dotted': ':',
    'dashdot': '-.'
}
