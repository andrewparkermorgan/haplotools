#! /usr/bin/env python

from __future__ import print_function

import os
import sys
import csv
import collections

from genome import Genome
genomes = {}

already = collections.defaultdict(str)

with open("/db/populations/CC/all.haps.csv", "r") as infile:
	lines = csv.DictReader(infile)
	## line format: line,id,name,phase,chromosome,start,end,strain,timestamp
	for line in lines:
		iid = line["line"]
		if iid in already:
			if already[iid] != line["id"]:
				continue
		else:
			already[iid] = line["id"]
		phase = 0 if line["phase"] == "top" else 1
		if not iid in genomes:
			genomes[iid] = Genome(iid, 0)
		genomes[iid].add_segment(line["chromosome"], int(line["start"]), int(line["end"]), phase, line["strain"], check = False)
		#print("{} -- {}:{}-{} [{}:{}]".format(line["line"],line["chromosome"], line["start"], line["end"], line["phase"], line["strain"]))

for iid in genomes.keys():
	print(iid, " ".join(genomes[iid].get_genotype_at("chr2", 83900000)))