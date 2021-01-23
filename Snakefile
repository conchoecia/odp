"""
This script performs mutual-best protein diamond BLAST searches,
 then makes synteny plots of those results.
"""
import matplotlib
import pandas as pd
import numpy as np
import sys

configfile: "config.yaml"

#make fake breaks
for this_axis in ["xaxisspecies", "yaxisspecies"]:
    for this_one in config[this_axis]:
        if "breaks" not in config[this_axis][this_one]:
            config[this_axis][this_one]["breaks"] = []


#make a function to add entries into y, if the entry is in x but not yet in y
for thisx in config["xaxisspecies"]:
    if thisx not in config["yaxisspecies"]:
        config["yaxisspecies"][thisx] = config["xaxisspecies"][thisx]

def expand_avoid_matching_x_and_y(filestring, xsamples, ysamples):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                outlist.append(filestring.format(xsamp, ysamp))
    return outlist

def expand_avoid_matching_x_and_y_third(filestring, xsamples, ysamples, third):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                for t in third:
                    outlist.append(filestring.format(xsamp, ysamp, t))
    return outlist


rule all:
    input:
        #expand("synteny_analysis/blastp_results/xtoybest/{xsample}_against_{ysample}.blastp",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"]),
        #expand("synteny_analysis/blastp_results/ytoxbest/{ysample}_against_{xsample}.blastp",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"])
        #expand("synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{ysample}_recip.blastp",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"]),
        #expand("synteny_analysis/genome_coords/x_genome_coords/{xsample}_genomecoords.txt",
        #       xsample = config["xaxisspecies"]),
        #expand("synteny_analysis/genome_coords/y_genome_coords/{ysample}_genomecoords.txt",
        #       ysample = config["yaxisspecies"])
        expand_avoid_matching_x_and_y("synteny_analysis/plots/synteny_uncolored/{}_and_{}_synteny.pdf",
                config["xaxisspecies"], config["yaxisspecies"]),
        #expand("synteny_analysis/plots/synteny_uncolored/{xsample}_and_{ysample}_synteny.pdf",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"]),
        expand_avoid_matching_x_and_y("synteny_analysis/dvalue_table/{}_and_{}_info.tsv",
                config["xaxisspecies"],  config["yaxisspecies"]),
        #expand("synteny_analysis/dvalue_table/{xsample}_and_{ysample}_info.tsv",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"]),

        # working on plotting for color
        #expand("synteny_analysis/prot_to_color/{colorby}_prottocolor.tsv",
        #       colorby = [x for x in config["xaxisspecies"] if "chrom_to_color" in config["xaxisspecies"][x]])
        expand_avoid_matching_x_and_y_third(
            "synteny_analysis/plots/synteny_colored_by/{}_and_{}_coloredby_{}_synteny.pdf",
            config["xaxisspecies"], config["yaxisspecies"],
            [x for x in config["xaxisspecies"] if "chrom_to_color" in config["xaxisspecies"][x]]
            ),
        #expand("synteny_analysis/plots/synteny_colored_by/{xsample}_and_{ysample}_coloredby_{colorby}_synteny.pdf",
        #       xsample = config["xaxisspecies"], ysample = config["yaxisspecies"],
        #       colorby = [x for x in config["xaxisspecies"] if "chrom_to_color" in config["xaxisspecies"][x]]),
        expand_avoid_matching_x_and_y_third(
            "synteny_analysis/plots/synteny_colored_by_no_missing/{}_and_{}_coloredby_{}_synteny.pdf",
            config["xaxisspecies"], config["yaxisspecies"],
            [x for x in config["xaxisspecies"] if "chrom_to_color" in config["xaxisspecies"][x]]
            ),
        #expand("synteny_analysis/plots/synteny_colored_by_no_missing/{xsample}_and_{ysample}_coloredby_{colorby}_synteny.pdf",
        #       xsample = config["xaxisspecies"], ysample = config["yaxisspecies"],
        #       colorby = [x for x in config["xaxisspecies"] if "chrom_to_color" in config["xaxisspecies"][x]]),
        # make protein similarity plots
        expand_avoid_matching_x_and_y(
            "synteny_analysis/plots/sample_similarity/{}_and_{}_peridentity_length.pdf",
             config["xaxisspecies"], config["yaxisspecies"]),
        #expand("synteny_analysis/plots/sample_similarity/{xsample}_and_{ysample}_peridentity_length.pdf",
        #        xsample = config["xaxisspecies"], ysample = config["yaxisspecies"]),

rule make_diamonddb_x:
    input:
        prots = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["proteins"]
    output:
        pep = "synteny_analysis/db/xaxis/{xsample}_prots.pep",
        dmnd = "synteny_analysis/db/xaxis/dmnd/{xsample}_prots.dmnd"
    threads: workflow.cores - 1
    shell:
        """
        ln -s {input.prots} {output.pep}
        diamond makedb --in {output.pep} --db {output.dmnd}
        """

rule make_diamonddb_y:
    input:
        prots = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["proteins"]
    output:
        pep = "synteny_analysis/db/yaxis/{ysample}_prots.pep",
        dmnd = "synteny_analysis/db/yaxis/dmnd/{ysample}_prots.dmnd"
    threads: workflow.cores - 1
    shell:
        """
        ln -s {input.prots} {output.pep}
        diamond makedb --in {output.pep} --db {output.dmnd}
        """

rule diamond_blast_x_to_y:
    input:
        xpep = "synteny_analysis/db/xaxis/{xsample}_prots.pep",
        ydmnd = "synteny_analysis/db/yaxis/dmnd/{ysample}_prots.dmnd",
    output:
        blastp = "synteny_analysis/blastp_results/xtoy/{xsample}_against_{ysample}.blastp",
    threads: workflow.cores - 1
    shell:
        """
        diamond blastp --query {input.xpep} --db {input.ydmnd} \
          --threads {threads} --evalue 1E-5 --outfmt 6 --out {output.blastp}
        """

rule diamond_blast_y_to_x:
    input:
        ypep = "synteny_analysis/db/yaxis/{ysample}_prots.pep",
        xdmnd = "synteny_analysis/db/xaxis/dmnd/{xsample}_prots.dmnd",
    output:
        blastp = "synteny_analysis/blastp_results/ytox/{ysample}_against_{xsample}.blastp",
    threads: workflow.cores - 1
    shell:
        """
        diamond blastp --query {input.ypep} --db {input.xdmnd} \
          --threads {threads} --evalue 1E-5 --outfmt 6 --out {output.blastp}
        """

rule absolute_best_from_blast_x_to_y:
    input:
        blastp = "synteny_analysis/blastp_results/xtoy/{xsample}_against_{ysample}.blastp",
    output:
        blastp = "synteny_analysis/blastp_results/xtoybest/{xsample}_against_{ysample}.blastp",
    threads: 1
    shell:
        """
        awk 'BEGIN{{former = ""}} {{if ($1 != former){{print($0)}}; former=$1}}' {input.blastp} > {output.blastp}
        """

rule absolute_best_from_blast_y_to_x:
    input:
        blastp = "synteny_analysis/blastp_results/ytox/{ysample}_against_{xsample}.blastp",
    output:
        blastp = "synteny_analysis/blastp_results/ytoxbest/{ysample}_against_{xsample}.blastp",
    threads: 1
    shell:
        """
        awk 'BEGIN{{former = ""}} {{if ($1 != former){{print($0)}}; former=$1}}' {input.blastp} > {output.blastp}
        """

rule reciprocal_best_hits:
    """
    finds the reciprocal best hits.
    reports it in the form of the blastp results from x -> y search
    """
    input:
        ytoxblastp = "synteny_analysis/blastp_results/ytoxbest/{ysample}_against_{xsample}.blastp",
        xtoyblastp = "synteny_analysis/blastp_results/xtoybest/{xsample}_against_{ysample}.blastp",
    output:
        xtoyblastp = "synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{ysample}_recip.blastp",
    threads: 1
    run:
        pairs = set()
        #first look in y_to_x
        with open(input.ytoxblastp, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split("\t")
                    add_this = (splitd[1], splitd[0])
                    if add_this in pairs:
                        raise IOError("This set was already found")
                    else:
                        pairs.add(add_this)
        # now go through x_to_y and filter
        out_handle = open(output.xtoyblastp, "w")
        with open(input.xtoyblastp, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split("\t")
                    check_this = (splitd[0], splitd[1])
                    if check_this in pairs:
                        print(line, file=out_handle)
        out_handle.close()

rule make_recip_table_for_marginplot:
    input:
        xtoyblastp = "synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{ysample}_recip.blastp",
    output:
        table = "synteny_analysis/blastp_results/reciprocal_best/table_for_marginplot/{xsample}_and_{ysample}_peridentity_length.tsv",
    threads: 1
    shell:
        """
        echo "" | awk '{{printf("PercentIdentity\\tMatchLength\\n")}}' > {output.table}
        cat {input.xtoyblastp} | cut -f3,4 >> {output.table}
        """

rule make_identity_marginplot:
    input:
        table = "synteny_analysis/blastp_results/reciprocal_best/table_for_marginplot/{xsample}_and_{ysample}_peridentity_length.tsv",
    output:
        plot = "synteny_analysis/plots/sample_similarity/{xsample}_and_{ysample}_peridentity_length.pdf"
    params:
        stem = "synteny_analysis/plots/sample_similarity/{xsample}_and_{ysample}_peridentity_length",
        xsample = lambda wildcards: wildcards.xsample,
        ysample = lambda wildcards: wildcards.ysample
    threads: 1
    shell:
        """
        pauvre custommargin -i {input.table} --fileform pdf \
          --xcol MatchLength --ycol PercentIdentity \
          -o {params.stem} --no_timestamp \
          --plot_min_y 0 --plot_max_y 100 \
          --plot_min_x 0 --plot_max_x 2000 \
          -t "Protein Identity of {params.xsample} and {params.ysample}"
        """

rule get_genome_coords_x:
    input:
        genome = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["genome"]
    output:
        coords = "synteny_analysis/genome_coords/x_genome_coords/{xsample}_genomecoords.txt",
    threads:
        1
    params:
        minsize = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["minscafsize"]
    shell:
        """
        bioawk -cfastx '{{ if (length($seq) >= {params.minsize}) {{ \
                           print($name, length($seq), sum)  }} \
                        }}' {input.genome} | \
          sort -k2 -nr | \
          awk '{{sum = sum + $2; print($1, $2, sum) }}' > {output.coords}
        """

rule get_genome_coords_y:
    input:
        genome = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["genome"]
    output:
        coords = "synteny_analysis/genome_coords/y_genome_coords/{ysample}_genomecoords.txt",
    threads:
        1
    params:
        minsize = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["minscafsize"]
    shell:
        """
        bioawk -cfastx '{{ if (length($seq) >= {params.minsize}) {{ \
                           print($name, length($seq), sum)  }} \
                        }}' {input.genome} | \
          sort -k2 -nr | \
          awk '{{sum = sum + $2; print($1, $2, sum) }}' > {output.coords}
        """
def blast_plot_order_helper(coords, sample, breaks, xory, xprottoloc, yprottoloc, recip, xorder):
    """
    This uses the reciprocal blast results to come up with the sort order
     for the y-axis scaffolds. Returns a list of the plot order.

    This code is all duplicated from the synteny plot function.
     Could be programmed in a better way to avoid redundancy, but this just fits
     the edge case where the y-axis has to be arranged based on the blast results.
    """
    # now make a lookup table of where the prots are.
    #  Use the x_offset and y_offset to recalculate where the plotting
    #  value is
    x_prot_to_loc = {}
    x_prot_to_scaf = {}
    with open(xprottoloc, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                prot = splitd[0]
                scaf = splitd[1]
                pos = int(splitd[2])
                x_prot_to_loc[prot] = pos
                x_prot_to_scaf[prot] = scaf
    y_prot_to_loc = {}
    y_prot_to_scaf = {}
    with open(yprottoloc, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                prot = splitd[0]
                scaf = splitd[1]
                pos = int(splitd[2])
                y_prot_to_loc[prot] = pos
                y_prot_to_scaf[prot] = scaf

    # now open the blast results and translate the pairs
    #  into plotting positions
    df = pd.read_csv(recip, header=None, sep = "\t")
    df.columns = ["xgene", "ygene", "pident", "length",
                  "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    df = df[["xgene", "ygene", "bitscore", "evalue"]]
    #print(x_prot_to_loc)
    df["xpos"] = df["xgene"].map(x_prot_to_loc)
    df["ypos"] = df["ygene"].map(y_prot_to_loc)

    df["yscaf"] = df["ygene"].map(y_prot_to_scaf)
    df["xscaf"] = df["xgene"].map(x_prot_to_scaf)
    df = df.dropna()
    df = df.sort_values(by=['xpos'])
    df = df.loc[df["evalue"] <= float("1E-20"), ]
    df = df.dropna()

    grouped_df = df.groupby(["yscaf"])
    for key, item in grouped_df:
        max_item = grouped_df.get_group(key)['xscaf'].value_counts().idxmax()
        all_other_things = [x for x in grouped_df.get_group(key)['xscaf'].unique() if x != max_item]
        for thisthing in all_other_things:
            df = df.loc[~( (df["yscaf"] == key) & (df["xscaf"] == thisthing)), ]
    # now sort based on the xscafs and the xpos
    sorterIndex = dict(zip(xorder, range(len(xorder))))
    df.sort_values(['yscaf', 'ypos'],
        ascending = [True, True], inplace = True)
    df.reset_index(drop=True, inplace = True)
    df = df.drop_duplicates(subset=['yscaf'])
    df['x_Rank'] = df['xscaf'].map(sorterIndex)
    df.sort_values(['x_Rank', 'xpos'],
        ascending = [True, True], inplace = True)
    df = df.dropna()
    df.reset_index(drop=True, inplace = True)
    print(df)
    #print(list(df.yscaf))
    return(list(df.yscaf))

def parse_coords(coords, sample, breaks, xory, xprottoloc=None, yprottoloc=None, recip=None, xorder=None):
    """
    This parses the coordinates and returns a
      - coord-to-offset dict (I don't remember what this is for),
      - a list of locations to plot lines (These are the scaf/chrom divisions)
      - the max value for that axis
      - the BOS positions for that to plot (dotted, translucent line)
      - the tick labels
      - the tick positions
      - the yorder or xorder
    """
    offset = {}
    max_coord = 0
    lines_at = []
    df = pd.read_csv(coords, header = None, sep = " ")
    df.columns = ["scaf", "scaflen", "cumsum"]
    # now figure out if we need to sort or not
    if xory == "x":
        if sample not in config["xaxisspecies"]:
            raise IOError("Can't find this xspecies")
        else:
            if "plotorder" in config["xaxisspecies"][sample]:
                plotorder = config["xaxisspecies"][sample]["plotorder"]
                print("plot order is: {}".format(plotorder))
            else:
                plotorder = None
    elif xory == "y":
        #print("we're in y")
        if sample not in config["yaxisspecies"]:
            raise IOError("Can't find this yspecies")
        else:
            #print("we're in the else of y")
            if "plotorder" in config["yaxisspecies"][sample]:
                print("we're in the plotorder of y")
                plotorder = config["yaxisspecies"][sample]["plotorder"]
            if "sort_by_x_coord_blast" in config["yaxisspecies"][sample]:
                #print("we're in the sort_by_x_coord_blast of y")
                if config["yaxisspecies"][sample]["sort_by_x_coord_blast"]:
                    #print("we're in the sort_by_x_coord_blast of y True")
                    # we need to set up the sort order based on the occurrence in the blast results
                    plotorder = blast_plot_order_helper(coords, sample, breaks,
                                                        xory, xprottoloc,
                                                        yprottoloc, recip, xorder)
                    #print("after plot order in sort_by_x_coord_blast")
            else:
                plotorder = None
    else:
        raise IOError("Don't know what this is")
    # now we have determined if we need to sort
    if plotorder != None:
        #print(" - using custom plot order: ", plotorder)
        sortdict = {key: val for key, val in zip(plotorder, range(len(plotorder)))}
        df['rank'] = df['scaf'].map(sortdict)
        df.sort_values(by = 'rank' ,inplace=True)
    df = df.dropna()
    df.reset_index(drop=True, inplace = True)
    df["cumsum"] = df["scaflen"].cumsum()
    df["cumsum"] = df["cumsum"] - df["scaflen"]
    print(df)
    for i, row in df.iterrows():
        offset[row["scaf"]] = row["cumsum"]
        if i > 0:
            lines_at.append(row["cumsum"])
    max_coord = list(df["scaflen"].cumsum())[-1]

    # now get the break of synteny positions
    BOS_here = []
    for pair in breaks:
        scaf = pair[0]
        pos  = pair[1]
        BOS_here.append( offset[scaf] + pos )

    #tick labels
    tick_labels = list(df["scaf"])
    tick_pos    = list(df["cumsum"] + (df["scaflen"]/2))
    return (offset, lines_at, max_coord, BOS_here, tick_labels, tick_pos, list(df["scaf"]))

def calc_D_for_y_and_x(df):
    """
    This calculates D for both the x and y axes.
    Defined in the 2020 vertebrate synteny paper.
    """
    df = df.dropna()
    df = df.sort_values(by=['xpos'])
    df.reset_index(drop=True, inplace = True)

    # this just calculates Dx
    df2 = pd.get_dummies(df["xgene_on_which_y_scaf"])
    # tried to weight the values but didn't work well
    #df2 = df.pivot_table(values='bitscore', index=df.index, columns='xgene_on_which_y_scaf', fill_value=0)
    #df2 = df2/df["bitscore"].max()
    #print(df2)
    #sys.exit()
    df2_xiL = df2.apply(lambda x: x.rolling(20).mean(), axis = 0)
    df2_xiR = df2.apply(lambda x: x.iloc[::-1].rolling(20).mean(), axis = 0).iloc[::-1]
    df2_xiR = df2_xiR.set_index(df2_xiR.index - 1)
    df2_xiR = df2_xiR.iloc[1:]
    subtractdf = df2_xiR.fillna(0) - df2_xiL.fillna(0)
    D = subtractdf.apply(lambda x: np.sqrt(np.square(x).sum()), axis = 1)
    df["Dx"] = D
    df["Dx_barwidth"] = [(list(df["xpos"])[i+1] - list(df["xpos"])[i])*0.8 for i in range(len(list(df["Dx"]))-1)] + [1]

    # this just calculates Dy
    df = df.sort_values(by=['ypos'])
    df.reset_index(drop=True, inplace = True)
    df2 = pd.get_dummies(df["ygene_on_which_x_scaf"])
    df2_yiL = df2.apply(lambda x: x.rolling(20).mean(), axis = 0)
    df2_yiR = df2.apply(lambda x: x.iloc[::-1].rolling(20).mean(), axis = 0).iloc[::-1]
    df2_yiR = df2_yiR.set_index(df2_yiR.index - 1)
    df2_yiR = df2_yiR.iloc[1:]
    subtractdf_y = df2_yiR.fillna(0) - df2_yiL.fillna(0)
    D = subtractdf_y.apply(lambda x: np.sqrt(np.square(x).sum()), axis = 1)
    df["Dy"] = D
    df["Dy_barwidth"] = [(list(df["ypos"])[i+1] - list(df["ypos"])[i])*0.8 for i in range(len(list(df["Dy"]))-1)] + [1]

    df = df.sort_values(by=['xpos'])
    df.reset_index(drop=True, inplace = True)
    return df

def synteny_plot(ycoords, xcoords, xprottoloc, yprottoloc, xsample, ysample,
                 xbreaks, ybreaks, recip, synplot, outtable, plotorder_file, prot_to_color,
                 dropmissing, plot_y_lines = False):
        import pandas as pd
        import seaborn as sns; sns.set()
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        from matplotlib.ticker import StrMethodFormatter, NullFormatter
        import numpy as np
        # set seaborn stuff
        #sns.set(rc={'text.usetex' : True})
        sns.set_style("ticks", {'font.family': ['sans-serif'],
                                    'font.sans-serif': ['Helvetica'],
                                    'grid.color': '.95'})
        # Preserve the vertical order of embedded images:
        matplotlib.rcParams['image.composite_image'] = False
        # text as font in pdf
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        # parse the xbreaks and the ybreaks.
        # These are where there are low-opacity,
        #  dotted lines to show breaks in synteny
        # example entries are like this:
          #"CM022887.1_chromosome_2:4392062"
          #"CM022887.1_chromosome_2:8784125"

        xbreaks_list = []
        for thisx in xbreaks:
            scaf = thisx.split(":")[0]
            pos  = int(thisx.split(":")[1])
            xbreaks_list.append([scaf, pos])
        ybreaks_list = []
        for thisy in ybreaks:
            scaf = thisy.split(":")[0]
            pos  = int(thisy.split(":")[1])
            ybreaks_list.append([scaf, pos])


        # first make a lookup table of how to calculate the
        #  x and y coords. This lookup is just the amount of
        # bp to add to the value when plotting. We pass the xprot_to_loc,
        #  xprot_to_scaf in case we need to sort everything based on order of
        #  occurrence on the scaffolds
        x_offset, vertical_lines_at, xmax, xbreaks, xticklabel, xtickpos, xorder = parse_coords(
            xcoords, xsample, xbreaks_list, "x")
        print("found {} x chromosomes".format(len(x_offset)))

        y_offset, horizontal_lines_at, ymax, ybreaks, yticklabel, ytickpos, yorder = parse_coords(
            ycoords, ysample, ybreaks_list, "y",
            xprottoloc, yprottoloc, recip, xticklabel)
        print("found {} y chromosomes".format(len(y_offset)))

        # now save the plot order to a file
        with open(plotorder_file, "w") as f:
            print("xplotorder:", file=f)
            for entry in xorder:
                print("  - {}".format(entry), file=f)
            print("yplotorder:", file=f)
            for entry in yorder:
                print("  - {}".format(entry), file=f)

        # now make a lookup table of where the prots are.
        #  Use the x_offset and y_offset to recalculate where the plotting
        #  value is
        x_prot_to_loc = {}
        x_prot_to_scaf = {}
        with open(xprottoloc, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split()
                    prot = splitd[0]
                    scaf = splitd[1]
                    if scaf in x_offset:
                        pos = int(splitd[2]) + x_offset[scaf]
                        x_prot_to_loc[prot] = pos
                        x_prot_to_scaf[prot] = scaf
        y_prot_to_loc = {}
        y_prot_to_scaf = {}
        with open(yprottoloc, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split()
                    prot = splitd[0]
                    scaf = splitd[1]
                    if scaf in y_offset:
                        pos = int(splitd[2]) + y_offset[scaf]
                        y_prot_to_loc[prot] = pos
                        y_prot_to_scaf[prot] = scaf

        # now open the blast results and translate the pairs
        #  into plotting positions
        df = pd.read_csv(recip, header=None, sep = "\t")
        df.columns = ["xgene", "ygene", "pident", "length",
                      "mismatch", "gapopen", "qstart", "qend",
                      "sstart", "send", "evalue", "bitscore"]
        df = df[["xgene", "ygene", "bitscore", "evalue"]]
        #print(x_prot_to_loc)
        df["xpos"] = df["xgene"].map(x_prot_to_loc)
        df["ypos"] = df["ygene"].map(y_prot_to_loc)
        print(df)

        df["xgene_on_which_y_scaf"] = df["ygene"].map(y_prot_to_scaf)
        df["ygene_on_which_x_scaf"] = df["xgene"].map(x_prot_to_scaf)
        df = df.dropna()
        df = df.sort_values(by=['xpos'])
        df = df.loc[df["evalue"] <= float("1E-20"), ]
        df.reset_index(drop=True, inplace = True)

        # Now calculate D and Dw for both X and Y axes
        df = calc_D_for_y_and_x(df)
        print(df)

        xgene = df["xgene"]
        x = df["xpos"]
        y = df["ypos"]
        bitscore = df["bitscore"]

        print("found {} points to plot".format(len(df["xpos"])))
        print("max bitscore: ", max(bitscore))
        bitscore_adjusted = [(x/max(bitscore))*60 for x in bitscore]
        colors = [(0, 0, 1.0, min( 1.0, (x/max(bitscore))*(max(bitscore)/np.mean(bitscore)))) for x in bitscore]
        drops = set()
        if prot_to_color:
            for i in range(len(xgene)):
                alpha = colors[i][3]
                try:
                    newcolor = list(matplotlib.colors.to_rgba(prot_to_color[xgene[i]]))
                except:
                    print("couldn't find a color for: {}".format(xgene[i]))
                    newcolor = [0,0,0,1]
                    drops.add(i)
                newcolor[3] = alpha
                colors[i] = newcolor

        if dropmissing:
            vararray = [xgene, x, y, bitscore, bitscore_adjusted, colors]
            for j in range(len(vararray)):
                vararray[j] = [vararray[j][i] for i in range(len(vararray[j])) if i not in drops]
            xgene = vararray[0]
            x     = vararray[1]
            y     = vararray[2]
            bitscore = vararray[3]
            bitscore_adjusted = vararray[4]
            colors = vararray[5]

        # save the output
        if outtable:
            df.to_csv(outtable, sep="\t")

        # now make a scatter plot
        figWidth = 8
        figHeight = 8
        plt.figure(figsize=(figWidth,figHeight))
        #set the panel dimensions
        panelWidth = 4
        panelHeight = 4
        dpanel_width = 0.25
        #find the margins to center the panel in figure
        leftMargin = (figWidth - panelWidth)/2
        bottomMargin = ((figHeight - panelHeight)/2)
        panel1 = plt.axes([leftMargin/figWidth, #left
                             bottomMargin/figHeight,    #bottom
                             panelWidth/figWidth,   #width
                             panelHeight/figHeight])     #height
        panelxd = plt.axes([leftMargin/figWidth, #left
                             (bottomMargin+panelHeight+0.1)/figHeight,    #bottom
                             panelWidth/figWidth,   #width
                             dpanel_width/figHeight])     #height
        panelyd = plt.axes([(leftMargin+panelWidth + 0.1)/figWidth, #left
                             bottomMargin/figHeight,    #bottom
                             dpanel_width/figWidth,   #width
                             panelHeight/figHeight])     #height
        panel1.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=True,
                            left=False, labelleft=True,
                            right=False, labelright=False,
                            top=False, labeltop=False)
        panelxd.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=False, labelright=False,
                            top=False, labeltop=False)
        panelyd.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=False, labelright=False,
                            top=False, labeltop=False)
        # set the panel linewidth thinner
        for this_panel in [panel1, panelxd, panelyd]:
            for axis in ['top','bottom','left','right']:
                this_panel.spines[axis].set_linewidth(0.5)
        # turn off the axis spines
        for this_panel in [panelxd, panelyd]:
            this_panel.spines['top'].set_visible(False)
            this_panel.spines['right'].set_visible(False)

        panel1.scatter(x, y, color = colors,
                       ec = None, s=6, linewidths = 0)
        # set mins and max
        panel1.set_xlim([0,xmax])
        panel1.set_ylim([0,ymax])
        # set x ticks
        panel1.set_xticks(xtickpos)
        panel1.set_xticklabels(xticklabel, rotation=90, fontsize=8)
        # set y ticks
        if not plot_y_lines:
            #there are inevitably going to be many scaffolds. We need to subset
            # get a list of evenly spaced indices
            numElems = 20
            arr = ytickpos
            idx = np.round(np.linspace(0, len(arr) - 1, numElems)).astype(int)
            newarr       = [arr[i] for i in idx]
            newarrlabels = [round(arr[i]/1000000, 1) for i in idx]
            newarrDy     = [i for i in idx]
            # turn on y-axis ticks
            panel1.tick_params(left=True)
            panel1.set_yticks(newarr)
            panel1.set_yticklabels(newarrlabels, fontsize=8)
            # turn on y-axis ticks on the Dy plot
            panelyd.tick_params(right=True, labelright=True)
            panelyd.set_yticks(newarr)
            panelyd.set_yticklabels(newarrDy, fontsize=8)
            panelyd.yaxis.set_label_position("right")
            panelyd.set_ylabel("number of scaffolds")
            #panel1.set_yticklabels(yticklabel, fontsize=8)
            panel1.set_ylabel(ysample + " Mb")
        else:
            panel1.set_yticks(ytickpos)
            panel1.set_yticklabels(yticklabel, fontsize=8)
            panel1.set_ylabel(ysample)
        #print(list(zip(df["xpos"], df["Dx"])))

        # set the x and y labels
        panel1.set_xlabel(xsample)

        panelxd.bar(x = df["xpos"], height=df["Dx"], width = df["Dx_barwidth"], lw=0, color="blue", zorder = 2)
        panelxd.set_xlim([0,xmax])
        panelxd.set_ylabel('Dx', fontsize=10)

        panelyd.barh(y = df["ypos"], width=df["Dy"], height = df["Dy_barwidth"], lw=0, color="blue", zorder = 2)
        panelyd.set_ylim([0,ymax])
        panelyd.set_xlabel('Dy', fontsize=10)

        for this_axis in [panel1, panelxd, panelyd]:
            this_axis.xaxis.get_offset_text().set_visible(False)
            this_axis.yaxis.get_offset_text().set_visible(False)

        #plot vertical lines
        for value in vertical_lines_at:
            panel1.axvline(x=value, color="black", lw=0.5)
        #plot horizontal lines
        if plot_y_lines:
            for value in horizontal_lines_at:
                panel1.axhline(y=value, color="black", lw=0.5)

        # plot vertical BOS
        for value in xbreaks:
            panel1.axvline(x=value, color=[0,0,0,0.25], lw=0.5, linestyle="dotted")
        # plot horizontal BOS
        for value in ybreaks:
            panel1.axhline(y=value, color=[0,0,0,0.25], lw=0.5, linestyle="dotted")
        plt.savefig(synplot)

"""
This makes the synteny plot without doing any special coloring of the dots
"""
rule plot_synteny:
    input:
        ycoords = "synteny_analysis/genome_coords/y_genome_coords/{ysample}_genomecoords.txt",
        xcoords = "synteny_analysis/genome_coords/x_genome_coords/{xsample}_genomecoords.txt",
        xprottoloc = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["prot_to_loc"],
        yprottoloc = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["prot_to_loc"],
        recip = "synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{ysample}_recip.blastp"
    output:
        synplot = "synteny_analysis/plots/synteny_uncolored/{xsample}_and_{ysample}_synteny.pdf",
        table = "synteny_analysis/dvalue_table/{xsample}_and_{ysample}_info.tsv",
        plot_order = "synteny_analysis/plot_order/basic/{xsample}_and_{ysample}_plotorder.tsv",
    threads:
        1
    params:
        xsample = lambda wildcards: wildcards.xsample,
        ysample = lambda wildcards: wildcards.ysample,
        xbreaks  = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["breaks"],
        ybreaks  = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["breaks"],
        keep_y   = lambda wildcards: False if "sort_by_x_coord_blast" in config["yaxisspecies"][wildcards.ysample] else True
    run:
        synteny_plot(input.ycoords, input.xcoords,
                     input.xprottoloc, input.yprottoloc,
                     params.xsample, params.ysample,
                     params.xbreaks, params.ybreaks,
                     input.recip, output.synplot, output.table,
                     output.plot_order, None,
                     False, params.keep_y)

rule xprot_to_color:
    """
    In addition to plotting the synteny without a color scheme,
      we also would like to plot by coloring with another species' color scheme
    The output is just:
    xsample_prot\thex_color
    """
    input:
        x_prot_to_loc = lambda wildcards: config["xaxisspecies"][wildcards.colorby]["prot_to_loc"]
    output:
        prot_to_color = "synteny_analysis/prot_to_color/{colorby}_prottocolor.tsv"
    params:
        colormap = lambda wildcards: config["xaxisspecies"][wildcards.colorby]["chrom_to_color"]
    run:
        print(params.colormap)
        # parse the printing information
        print_list = []
        for key in params.colormap:
            coord = key
            color = params.colormap[coord]
            scaf = coord.split(":")[0]
            pos_raw = coord.split(":")[1]
            if pos_raw == "all":
                pos_min = 1
                pos_max = 999999999
            else:
                pos_min = int(pos_raw.split("-")[0])
                pos_max = int(pos_raw.split("-")[1])
            print_list.append({"scaf": scaf, "pos_min": pos_min,
                               "pos_max": pos_max, "color": color})
        df = pd.DataFrame.from_dict(print_list)
        print(df)
        #df["colors_py"] =  df["color"].apply(matplotlib.colors.to_rgba)
        out_handle = open(output.prot_to_color, "w")
        with open(input.x_prot_to_loc, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split()
                    prot = splitd[0]
                    scaf = splitd[1]
                    if scaf in list(df["scaf"]):
                        pos  = int(splitd[2])
                        query = "scaf == '{}' & pos_min <= {} & pos_max >= {}".format(scaf, pos, pos)
                        color = df.query(query)["color"].values[0]
                        print("{}\t{}".format(prot, color), file=out_handle)
                    else:
                        color = "#000000"
                        print("{}\t{}".format(prot, color), file=out_handle)
        out_handle.close()

"""
This makes the plots where the proteins are colored by another protein set.

Jan 12 2021 I am not sure what the no missing plot version is.

For some reason I am getting an InputFunctionException
"""
rule plot_synteny_x_colored_by_x:
    input:
        ycoords = "synteny_analysis/genome_coords/y_genome_coords/{ysample}_genomecoords.txt",
        xcoords = "synteny_analysis/genome_coords/x_genome_coords/{xsample}_genomecoords.txt",
        xprottoloc = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["prot_to_loc"],
        yprottoloc = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["prot_to_loc"],
        recip = "synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{ysample}_recip.blastp",
        color_by_blast = "synteny_analysis/blastp_results/reciprocal_best/{xsample}_and_{colorby}_recip.blastp",
        color_by = "synteny_analysis/prot_to_color/{colorby}_prottocolor.tsv"
    output:
        synplot = "synteny_analysis/plots/synteny_colored_by/{xsample}_and_{ysample}_coloredby_{colorby}_synteny.pdf",
        nodots = "synteny_analysis/plots/synteny_colored_by_no_missing/{xsample}_and_{ysample}_coloredby_{colorby}_synteny.pdf",
        plot_order = "synteny_analysis/plot_order/colored_by/{xsample}_and_{ysample}_coloredby_{colorby}_plotorder.tsv"
    threads:
        1
    params:
        xsample = lambda wildcards: wildcards.xsample,
        ysample = lambda wildcards: wildcards.ysample,
        color_sample = lambda wildcards: wildcards.colorby,
        xbreaks  = lambda wildcards: config["xaxisspecies"][wildcards.xsample]["breaks"],
        ybreaks  = lambda wildcards: config["yaxisspecies"][wildcards.ysample]["breaks"],
        keep_y   = lambda wildcards: False if "sort_by_x_coord_blast" in config["yaxisspecies"][wildcards.ysample] else True
    run:
        #first figure out what the color should be for each protein pair
        colorby_prot_to_color = {}
        with open(input.color_by, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    splitd = line.split("\t")
                    xprot = splitd[0]
                    colorbyprot = splitd[1]
                    colorby_prot_to_color[xprot] = colorbyprot
        # now use the blast results to lookup the prot color
        if params.xsample == params.color_sample:
            print("x is the same as color")
            prot_to_color = colorby_prot_to_color
        else:
            prot_to_color = {}
            try_counter =0
            fail_counter = 0
            with open(input.color_by_blast, "r") as f:
                for line in f:
                    line = line.strip()
                    if line:
                        splitd = line.split("\t")
                        try:
                            prot_to_color[splitd[0]] = colorby_prot_to_color[splitd[1]]
                            try_counter += 1
                        except:
                            prot_to_color[splitd[0]] = "#000000"
                            fail_counter += 1
            print(try_counter, fail_counter)
        # now that we have the colors, plot
        synteny_plot(input.ycoords, input.xcoords,
                     input.xprottoloc, input.yprottoloc,
                     params.xsample, params.ysample,
                     params.xbreaks, params.ybreaks,
                     input.recip, output.synplot, None,
                     output.plot_order, prot_to_color,
                     False, params.keep_y)
        # plot again, but without the black (missing) points
        synteny_plot(input.ycoords, input.xcoords,
                     input.xprottoloc, input.yprottoloc,
                     params.xsample, params.ysample,
                     params.xbreaks, params.ybreaks,
                     input.recip, output.nodots, None,
                     output.plot_order, prot_to_color,
                     True, params.keep_y)
