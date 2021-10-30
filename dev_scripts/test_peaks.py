import pandas as pd
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter


def flatten(list_of_lists):
   return [item for sublist in list_of_lists for item in sublist]
    
df = pd.read_csv("/Users/darrin/Downloads/caps_and_rhopilema_info.tsv",
                    sep = "\t", index_col = 0)
sort_order = {"x": {"pos": "xpos", "chrom": "ygene_on_which_x_scaf", "D": "Dx"}, 
              "y": {"pos": "ypos", "chrom": "xgene_on_which_y_scaf", "D": "Dy"}}
sort_direction = "x"

df = df.sort_values(by=[sort_order[sort_direction]["pos"]])
unique_chroms = []
chrom_breakpoints = []
for index, row in df.iterrows():
    thischrom = row[sort_order[sort_direction]["chrom"]]
    thispos  = row[sort_order[sort_direction]["pos"]]
    if thischrom not in unique_chroms:
        unique_chroms.append(thischrom)
        chrom_breakpoints.append(thispos)
chrom_breakpoints = chrom_breakpoints[1::]

df_lists = []

all_ranges = set()
for thiscol in ["deltMA","deltD", "Dx2"]:
    for thischrom in unique_chroms:
        subdf = df.loc[df[sort_order[sort_direction]["chrom"]] == thischrom, ].copy()
    
        subdf["Dx2"] = subdf["Dx"]
        subdf['Dx2'] = np.where((subdf["Dx"] < subdf["Dx"].median()),np.NaN,subdf["Dx2"])
        subdf["MA"] = subdf["Dx2"].rolling(window=3, center=True).mean()
        subdf["MA2"] = subdf["Dx2"].rolling(window=19, center=True).mean()   
    
        subdf["deltMA"] = subdf["MA"].diff() / subdf["MA"].index.to_series().diff()
        subdf['deltMA'] = np.where((subdf["MA"] < subdf["MA"].median()),np.NaN,subdf["deltMA"])
        
        subdf["deltD"] = subdf["Dx2"].diff() / subdf["Dx2"].index.to_series().diff()
        subdf['deltD'] = np.where((subdf.Dx2 < subdf["Dx2"].median()),np.NaN,subdf.deltD)
    
        idxmaxes = set()
        print(thischrom)
        # get the groups of consecutive values
        #for thiscol in ["deltMA","deltD", "Dx2"]:
        ind = list(subdf[~subdf[thiscol].isnull()].index)
        ranges =[]
        for k,g in groupby(enumerate(ind),lambda x:x[0]-x[1]):
            group = (map(itemgetter(1),g))
            group = list(map(int,group))
            ranges.append((group[0],group[-1]))   
    
        # now get the peak from this set of values
        print(thischrom, thiscol, ranges)
        if len(ranges) > 0:
            for this_range in ranges:
                if this_range[0] != this_range[-1]:
                    #this_range = [x for x in range(this_range[0], this_range[1]+1)]
                    this_range = list(this_range)
                    print("    {}".format(this_range))
                    which_d_col = sort_order[sort_direction]["D"]
                    temp = subdf.loc[this_range[0]:this_range[-1]][which_d_col].idxmax()
                    idxmaxes.add(temp)
    
        print(thischrom, idxmaxes)
        keep_idx_maxes = set()
        # picks the best in windows of 50 genes
        window = 20
        ignore_set = set()
        done = False
        consider_ranges = set()
        while not done:
            consider_ranges = set()
            for this_idx in idxmaxes:
                # get windows of ranges if they're not in the ignore set
                thistup = tuple([x for x in idxmaxes
                     if ((x > this_idx - window)
                         and (x < this_idx + window)
                         and (x not in ignore_set))])
                if len(thistup) > 0:
                    consider_ranges.add(thistup)
            
            consider_ranges = sorted(list(consider_ranges), key=len, reverse=True)
            if len(consider_ranges) > 0:
                thisrange = list(consider_ranges[0])
                if len(thisrange) == 1: 
                    done = True
                else:
                    submax = df.loc[thisrange, ]["Dx"].idxmax()
                    for thisid in thisrange:
                        if thisid != submax:
                            ignore_set.add(thisid)
            else:
                # there's nothing here
                done = True
        for entry in consider_ranges:
            all_ranges.add(entry)  
        print()

# flatten the results of what we got from the last analysis

idxmaxes = flatten(all_ranges) 
# now pick the best in windows of 10 genes after combining
window = 5
ignore_set = set()
done = False
consider_ranges = set()
while not done:
    consider_ranges = set()
    for this_idx in idxmaxes:
        # get windows of ranges if they're not in the ignore set
        thistup = tuple([x for x in idxmaxes
             if ((x > this_idx - window)
                 and (x < this_idx + window)
                 and (x not in ignore_set))])
        if len(thistup) > 0:
            consider_ranges.add(thistup)
    
    consider_ranges = sorted(list(consider_ranges), key=len, reverse=True)
    if len(consider_ranges) > 0:
        thisrange = list(consider_ranges[0])
        if len(thisrange) == 1: 
            done = True
        else:
            submax = df.loc[thisrange, ]["Dx"].idxmax()
            for thisid in thisrange:
                if thisid != submax:
                    ignore_set.add(thisid)
    else:
        # there's nothing here
        done = True

vert_lines = flatten(consider_ranges) 

plt.figure(figsize=(30, 6))
plt.plot(df["xpos"], df["Dx"])
#plt.plot(df["xpos"], df["Dx2"])
#plt.plot(df["xpos"], df["MA"], color = "red")
#plt.plot(df["xpos"], (df["delt"]*5) + df["MA"].mean(), color = "green")

# change the xlims
plot_max_x = max(df[sort_order[sort_direction]["pos"]])
plt.xlim([0, plot_max_x])

for thischrom in unique_chroms:
    subdf = df.loc[df[sort_order[sort_direction]["chrom"]] == thischrom, ]
    minx = min(subdf[sort_order[sort_direction]["pos"]])
    maxx = max(subdf[sort_order[sort_direction]["pos"]])
    print(minx, maxx)
    plt.axhline(y= subdf["Dx"].mean(),
                color='red',
                linestyle='-',
                xmin = minx/plot_max_x,
                xmax = maxx/plot_max_x)
                
    plt.axhline(y= subdf["Dx"].median(),
                color='green',
                linestyle='-',
                xmin = minx/plot_max_x,
                xmax = maxx/plot_max_x)

#plt.axhline(y= df["Dx"].median(),
#                color='red',
#                linestyle='-')

#for position in idxmaxes:
for position in vert_lines:
    xpos = df.loc[position]["xpos"]
    plt.axvline(x=xpos, color="black", ls = "--", lw = 0.75, alpha=0.5)

# plot the chrombreaks as vertical lines
for thischrom in chrom_breakpoints:
    plt.axvline(x=thischrom, color="black", lw = 2, alpha=1)


plt.show()