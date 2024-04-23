"""
These are some helper functions to make creating bokeh plots easier.
For example, to change colors of a plot, you need to convert the hex strings that include alpha to integers.
There are also multiple calls needed to remove the ticks from one plot.
"""
import numpy as np

def convert_hex_string_to_colorvalues(list_of_color_strings) -> np.ndarray:
    """
    Converts something like: ["#1F779A", "#9F0000"]
     to a numpy array like:
     [0x1F779AFF, 0x9F0000FF]
    The numpy array is what Bokeh expects for colors.
    """
    # first check every string starts with a #
    for color in list_of_color_strings:
        if color[0] != "#":
            raise ValueError(f"Color {color} does not start with a #")
    # then convert the strings to integers
    return np.array(
        [int(color[1:] + "FF", 16) for color in list_of_color_strings],
        dtype=np.uint32)

def remove_ticks(plot) -> None:
    """
    Removes the ticks from a Bokeh plot.
    """
    plot.xaxis.major_tick_line_color = None
    plot.xaxis.minor_tick_line_color = None
    plot.yaxis.major_tick_line_color = None
    plot.yaxis.minor_tick_line_color = None
    plot.xaxis.major_label_text_font_size = '0pt'
    plot.yaxis.major_label_text_font_size = '0pt'
