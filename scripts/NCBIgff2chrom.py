#!/usr/bin/env python3
"""
This program parses a NCBI GFF annotation and generates a .chrom file.
see https://github.com/conchoecia/odp for the specification.
"""
#NW_011887297.1	RefSeq	CDS	1566678	1566739	.	+	2	ID=cds-XP_004348322.1;Parent=rna-XM_004348272.1;Dbxref=GeneID:14898863,Genbank:XP_004348322.1;Name=XP_004348322.1;gbkey=CDS;locus_tag=CAOG_04494;product=hypothetical protein;protein_id=XP_004348322.1

import gzip
import sys

prots = {}

gzipped = False
for thisend in [".gz", ".gzip", ".GZ", ".GZIP", ".gzipped", ".GZIPPED"]:
    if sys.argv[1].endswith(thisend):
        gzipped = True

if gzipped:
    handle = gzip.open(sys.argv[1],'rt')
else:
    handle = open(sys.argv[1], "r")

for line in handle:
    line = line.strip()
    splitd=line.split("\t")
    if line and len(splitd) > 7 and splitd[2] == "CDS" and "protein_id=" in line:
        pid = [x for x in splitd[8].split(";") if x.startswith("protein_id=")][0].replace("protein_id=", "")
        scaf = splitd[0]
        strand = splitd[6]
        start = int(splitd[3])
        stop = int(splitd[4])
        if pid not in prots:
            prots[pid] = {"scaf": scaf, "strand": strand,
                          "start": start, "stop": stop}
        else:
            if start < prots[pid]["start"]:
                prots[pid]["start"] = start
            if stop > prots[pid]["stop"]:
                prots[pid]["stop"] = stop
handle.close()

for pid in prots:
    print("{}\t{}\t{}\t{}\t{}".format(
        pid, prots[pid]["scaf"],
        prots[pid]["strand"], prots[pid]["start"], prots[pid]["stop"]))
