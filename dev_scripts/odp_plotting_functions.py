"""
These are the plotting functions used by many ODP programs
"""

import matplotlib

def format_matplotlib():
    """format the fonts and print options for the plots"""
    font = {'family' : 'sans-serif',
            'sans-serif' : 'DejaVu Sans', # removing this after finding that many users don't have Helvetica installed. :( https://github.com/conchoecia/odp/issues/34
            'weight' : 'normal',
            'size'   : 12}

    matplotlib.rc('font', **font)

    grid = {"color": ".95", "linestyle": "-"}
    # grid style
    matplotlib.rc('grid', **grid)

    # Preserve the vertical order of embedded images:
    matplotlib.rcParams['image.composite_image'] = False
    # text as font in pdf
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
