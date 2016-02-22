#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import csv
import argparse

from hapfiles import HapfileParser

parser = argparse.ArgumentParser(description = "Convert csbio 'hapfile' format to a more sensible csv format.")
parser.add_argument(	"--recombs", action = "store_true",
						help = "get recombinations instead of haplotype blocks" )
parser.add_argument(	"infiles", nargs = "+",
						help = "paths to hapfiles to process" )
args = parser.parse_args()

# set up output stream
writer = csv.writer(sys.stdout, delimiter = "\t")

done_header = False
print("Parsing haplotypes for {} samples.".format(len(args.infiles)), file = sys.stderr)
for ff in args.infiles:

	if not done_header:
		#writer.writeheader()
		done_header = True

	sample = os.path.basename(ff).replace(".hap", "")
	sys.stderr.write("\t-- {}".format(sample))
	try:
		with open(ff, "r") as hapfile:
			parser = HapfileParser(hapfile)
			parser.parse()
			if not args.recombs:
				for b in parser.get_hapblocks():
					writer.writerow([ b.chrom, b.start, b.end, sample, b.founder, b.phase ])
			else:
				for b in parser.get_recombs():
					writer.writerow([ b.chrom, b.start, b.end, sample, b.attributes["from"], b.attributes["to"], b.phase ])
	except Exception as e:
			#sys.stderr.write(str(e) + "\n")
			sys.stderr.write( " --> FAILED")

	sys.stderr.write("\n")
