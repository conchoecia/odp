#!/usr/bin/env python

"""
This contains the helper scripts to annotate the UMAP plots with the blast results.
There is a module, because this allows us to import the functions into other programs for testing.
"""

from ast import literal_eval
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pandas as pd
import random

# imports form other scripts I have written:
from bokeh_helper import convert_hex_string_to_colorvalues, remove_ticks

def umapdf_one_species_one_query(UMAPdf, blastp, analysis, ALG, n, m, query, outputPDF, species = None):
    """
    This function takes a dataframe as input and makes a UMAP
        basedir + "/blast_filt/{query}/{sample}_results.filt.blastp",
        UMAPdf     = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        algrbhfile = config["ALG_rbh_file"],

    The UMAPdf is a dataframe with the positions, and the columns UMAP1, UMAP2, gene_group, color

    Output:
      - The output of this script is a UMAP pdf with the reembedding of the ALG, for the specified clade.
      - This also outputs a pandas dataframe that includes the closest blast results.
    """
    if species == None:
        raise IOError("Please use this function with a species. The inputs are different from the other function and this plots something in the context of one species. Exiting")
    if not outputPDF.endswith(".pdf"):
        raise IOError("The outputPDF must end with .pdf. Exiting")
    #         ┓   ┓•┓
    # ┏┳┓┏┓╋┏┓┃┏┓╋┃┓┣┓
    # ┛┗┗┗┻┗┣┛┗┗┛┗┗┗┗┛
    #       ┛
    UMAPempty   = False
    BLASTPempty = False
    if os.stat(UMAPdf).st_size == 0:
        UMAPempty = True
    if os.stat(blastp).st_size == 0:
        BLASTPempty = True
    if UMAPempty or BLASTPempty:
        # make a 2in x 2in plot telling the user the message
        fig = plt.figure(figsize=(2,2))
        if UMAPempty and BLASTPempty:
            plt.text(0.5, 0.5, "The input df and the blastp file were empty, so we couldn't plot anything", fontsize = 3)
        elif UMAPempty and not BLASTPempty:
            plt.text(0.5, 0.5, "The input df was empty, so we couldn't plot anything", fontsize = 3)
        elif not UMAPempty and BLASTPempty:
            plt.text(0.5, 0.5, "The blastp file was empty, so we couldn't plot anything", fontsize = 3)
        plt.text(0.5, 0.5, "The input df was empty, so we couldn't plot anything", fontsize = 3)
        # make another line telling the user which file, exactly, was empty
        plt.text(0.5, 0.4, f"The input file was {UMAPdf}", fontsize = 3)
        # turn off the axis ticks
        for ax in fig.axes:
            ax.set_xticks([])
            ax.set_yticks([])
        # turn off the axes
        plt.axis('off')
        try:
            # make the plot tight to not cut off the text
            plt.tight_layout()
        except:
            pass
        # save the figure
        plt.savefig(outputPDF)
    else:
        dot = 3
        bigger_dot = 4
        # The embedding in this case only has the UMAP1 and UMAP2 columns.
        # We need later to go through the blast results to figure out which dots to annotate
        df_embedding = pd.read_csv(UMAPdf, sep = "\t", index_col = 0)
        # add a column to the df_embedding that contains the BLASTP best hit, the distance of the hit to the ALG, and the color of the dot
        df_embedding["blastp_best_hits"] = [[] for i in range(df_embedding.shape[0])]
        df_embedding["blastp_distances"] = [[] for i in range(df_embedding.shape[0])]
        df_embedding["blastp_color"]     = ""
        print(df_embedding)
        # make a matplotlib plot of the UMAP with the df_embedding, and the color_dict from SplitLossColocTree as the legend
        # make a figure that is 5x5 inches
        fig, ax = plt.subplots(figsize=(2, 2))
        # scatter the UMAP1 and UMAP2 columns of the df_embedding
        ax.scatter(df_embedding["UMAP1"], df_embedding["UMAP2"], c = df_embedding["color"], s = dot)
        # now we read in the filtered blast results for this species
        blastp = pd.read_csv(blastp, sep = "\t", header = None)
        blastp.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore", "scaf_of_nearest_ALG", "nearest_ALG", "position_of_nearest_ALG", "dist_to_nearest_ALG"]
        blastp["closest_ALG"] = None
        blastp["closest_ALG_dist"] = None
        blastp["closest_ALG_position"] = None
        # get the blastp rows with nearest_ALG values in the UMAPdf.index
        for i, row in blastp.iterrows():
            thisquery = row["qseqid"]
            nearest_list = literal_eval(row["nearest_ALG"])
            for j in range(len(nearest_list)):
                if nearest_list[j] in df_embedding.index:
                    print(f"found a match for ALG {ALG}, blastp query {thisquery} {nearest_list[j]} in the UMAPdf.index")
                    nearest_pos = literal_eval(row["position_of_nearest_ALG"])[j]
                    nearest_dist = literal_eval(row["dist_to_nearest_ALG"])[j]
                    # We update the blastp dataframe with the nearest ALG, the distance, and the position. We use this for plotting dots later.
                    blastp.at[i, "closest_ALG"]          = nearest_list[j]
                    blastp.at[i, "closest_ALG_dist"]     = nearest_dist
                    blastp.at[i, "closest_ALG_position"] = nearest_pos
                    # We also update the dataframe with this info
                    df_embedding.at[nearest_list[j], "blastp_best_hits"].append(thisquery)
                    df_embedding.at[nearest_list[j], "blastp_distances"].append(nearest_dist)
                    break
        # remove rows where closest_ALG is None
        blastp = blastp[blastp["closest_ALG"].notnull()]

        # In many cases, the closest ALG is going to be the same for many of the blastp hits, so we should group these together
        gb = blastp.groupby("closest_ALG")
        # Go through each group, and color that dot in the UMAP
        legend_patches = []
        annotate_marker = '*'
        for thisALG, group in gb:
            # generate a random color for this dot
            thiscolor = "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            # add this color to the df_embedding
            df_embedding.at[thisALG, "blastp_color"] = thiscolor
            # label is all of the qseqids that are in this group
            label = ", ".join(group["qseqid"].values)
            legend_patches.append(  mlines.Line2D([], [], color=thiscolor, marker=annotate_marker, linestyle='None', markersize=dot*1.2, label=label))
            # x is the location of this ALG in the UMAPdf
            x = df_embedding.loc[thisALG, "UMAP1"]
            y = df_embedding.loc[thisALG, "UMAP2"]
            ax.scatter(x, y, c = thiscolor, s = dot, marker = annotate_marker, alpha = 1)

        # add the entries to the legend
        ax.legend(handles=legend_patches, loc="upper right", fontsize = 3)
        # Remove the frame from the legend
        ax.get_legend().set_frame_on(False)
        # turn off the axis ticks
        ax.set_xticks([])
        ax.set_yticks([])
        # set the title based on the input
        ax.set_title(f"UMAP of {analysis}, ALG {ALG},\ngenome {species} with {query} blast results,\nmin_dist {m}, n_neighbors {n}",
                   fontsize = 3)
        # save the figure
        plt.savefig(outputPDF, bbox_inches='tight')

        # Save the df_embedding with the blastp results, with the index saved
        df_embedding.to_csv(outputPDF.replace(".pdf", ".df"), sep = "\t", index = True)

import numpy as np
from bokeh import events
from bokeh.models import ColumnDataSource, CustomJS, TextInput, HoverTool, Div, Legend
from bokeh.plotting import figure, save, output_file
from bokeh.models.widgets import Button
# This sets out the layout to make a button to export the data.
from bokeh.layouts import column, row


def umapdf_reimbedding_bokeh_plot_one_species(blastdf, plot_title, output_html, scalar = 1.5):
    """
    This function takes one df from the reimbedding and makes a bokeh plot of the UMAP.
    The reasons for making these plots are that:
      - We will be able to hover over the dots and see the blast results
      - We will be able to use the lasso tool to select dots, and export a table of the selected dots.
    If we export these dots, we can use them to set intersection on other reimbeddings to see what points co-occur.

    Parameters:
      - blastdf: a pandas dataframe with the reimbedding and the blast results. The indices are the locus name.
                 These files are output by the snakemake rule `annotate_blastresults`
      - scalar:  a scalar to multiply the size of the plot. The default is 1.5, which makes the plot 600x600 pixels.
    """
    # check that blastdf is a file that exists
    if not os.path.exists(blastdf):
        raise IOError(f"The blastdf file {blastdf} does not exist. Exiting")
    # check that the type of plot_title is a string and not empty
    if (type(plot_title) != str) or (len(plot_title) == 0):
        raise IOError("The plot_title must be a non-empty string. Exiting")
    # first check that the output_html file path ends in .html
    if not output_html.endswith(".html"):
        raise IOError("The output_html file must end with .html. Exiting")
    # check that scalar is an int or a float
    if (type(scalar) != int) and (type(scalar) != float):
        raise IOError("The scalar must be an int or a float. Exiting")

    # Load the tsv. It has indices as the first column
    df = pd.read_csv(blastdf, sep="\t", index_col=0)

    # change the colors of the dots
    colorlist = []
    for i, thisrow in df.iterrows():
        # if blastp_color is NaN, then use the color column
        # Use a default blue if there is no color
        if pd.isna(thisrow['blastp_color']):
            if "#" not in thisrow['color']:
                colorlist.append("#1F779A")
            else:
                colorlist.append(thisrow['color'])
        else:
            if "#" not in thisrow['blastp_color']:
                colorlist.append("#9F0000")
            else:
                colorlist.append(thisrow['blastp_color'])
    colors = convert_hex_string_to_colorvalues(colorlist)

    s1 = ColumnDataSource(data=dict(locus       = list(df.index),
                                    UMAP1       = list(df["UMAP1"]),
                                    UMAP2       = list(df["UMAP2"]),
                                    blastpHits  = list(df['blastp_best_hits']),
                                    color       = colors,
                                    colorHex    = colorlist
                                    ))

    # Plot the data and save the html file
    p = figure(
        width  = int(400 * scalar),
        height = int(400 * scalar),
        tools  = ["lasso_select", "reset", "save"],
        title  = plot_title)
    p.circle('UMAP1', 'UMAP2',
             size=8, source=s1, alpha=0.4,
             line_color = None,
             fill_color="color",
             selection_color="firebrick")
    remove_ticks(p)

    # For the figure, add hover tools to see info about each point
    hover_tool = HoverTool(
        tooltips=[
            ("Locus",       "@locus"),
            ("UMAP1",       "@UMAP1"),
            ("UMAP2",       "@UMAP2"),
            ("Blastp best hits", "@blastpHits")
        ]
    )
    p.add_tools(hover_tool)

    # Define a button to change the save file name
    text_input = TextInput(value="selected-data.txt", title="Export Filename:")
    text_input.js_on_change("value", CustomJS(code="""
        console.log('text_input: value=' + this.value, this.toString());
    """))

    # Define a JavaScript callback to export the data
    # Create a Button widget for exporting data
    savebutton = Button(label="Export Data", button_type="success")
    savebutton.js_on_event(events.ButtonClick,
        CustomJS(args=dict(source_data=s1, filename=text_input),
        code="""
             var inds = source_data.selected.indices;
             var data = source_data.data;
             var out = "locus\\tUMAP1\\tUMAP2\\tcolor\\tblastp_best_hits\\n";
             for (var i = 0; i < inds.length; i++) {
                 out += data['locus'][inds[i]] + "\\t" + data['UMAP1'][inds[i]] + "\\t" + data['UMAP2'][inds[i]] + "\\t" + data['colorHex'][inds[i]] + "\\t" + data['blastpHits'][inds[i]] + "\\n";
             }
             var file = new Blob([out], {type: 'text/plain'});
             var elem = window.document.createElement('a');
             elem.href = window.URL.createObjectURL(file);
             elem.download = filename.value;
             document.body.appendChild(elem);
             elem.click();
             document.body.removeChild(elem);
             """
            )
        )

    # Underneath the button, add a label to give some instructions to the user.
    text_content="""
    <h1 style ="color:red;"> READ ME </h2>
    <p> The "Export Data" button above exports a dataframe of the lasso'd data.</p>
    <p>1. To use this plot, click on the lasso tool on the upper right of the plot.</p>
    <p>2. Then, type in a filename in the box above for your download.</p>
    <p>3. Click the "Export Data" button.</p>
    <p>\n</p>
    <p>This will export a file with a tab-separated list of the selected data points.</p>
    """
    text_box = Div(text=text_content, width=int(200 * scalar), height=int(200*scalar))

    # Show the plot and the button
    output_file(output_html)
    save(row(p, column(text_input,savebutton, text_box)))
