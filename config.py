"""
Visualization configuration for CIUT Riesgo project.
Contains color palettes, basemap settings, and line style definitions.
"""

import contextily as ctx

# Continuous color palettes
# Blue palette (5 sequential shades)
BLUE_PALETTE = [
    '#0d1b4d',  # Darkest blue
    '#1a2e5c',  # Dark blue
    '#2d4b73',  # Medium-dark blue
    '#4d7b99',  # Medium blue
    '#7ba3c2'   # Light blue
]

# Pink palette (5 sequential shades)
PINK_PALETTE = [
    '#4d0d1b',  # Darkest pink
    '#6a1a2e',  # Dark pink
    '#8a2d4b',  # Medium-dark pink
    '#a94d7b',  # Medium pink
    '#c27ba3'   # Light pink
]

# Binary palette: blue -> white -> pink
BINARY_PALETTE = ['#0d1b4d', '#ffffff', '#4d0d1b']

# Categorical color scheme (6 colors)
# Blue and pink from above + 4 colors that complement the scheme
CATEGORICAL_COLORS = [
    '#0d1b4d',  # Blue alta
    '#4d0d1b',  # Pink alta
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
