#! /usr/bin/env python

import os
import sys
import re
import csv
from collections import OrderedDict

from genome import Genome

class HapfileParser:

	def __init__(self, hapfile, sample = None):

		self._sample = sample
		self.file = hapfile
		self.hapblocks = None

	# does current line contain haplotype information or just display instructions?
	def _is_chromline(self, line):
		return ( re.match(r"^chr", line) is not None )

	# is this chromosome on maternal or paternal haplotype?
	# assume hapfile lists them as MATERNAL first
	def _get_phase(self, chrom, chroms_seen):
		if chrom in chroms_seen:
			return 1
		else:
			return 0

	# parse a line in the hap file into a list of blocks
	# note that the haplotype blocks are returned in the (native) csbio format:
	# |-----|
	#     |------|
	# whereas we want them to look like:
	# |---| |----|
	def _parse_chromline(self, line):
		blocks = []
		pieces = line.rstrip().split(",")
		chrom = pieces.pop(0)
		pieces.pop(0)
		founders = pieces[ ::3 ]
		starts = pieces[ 1::3 ]
		ends = pieces[ 2::3 ]
		new_starts, new_ends = [0]*len(starts), [0]*len(ends)
		new_starts[0] = starts[0]
		new_ends[ len(new_ends)-1 ] = ends[-1:][0]
		for i in range(1, len(starts)):
			new_starts[i] = ends[i-1]
			new_ends[i-1] = starts[i]
		for i in range(0, len(new_starts)):
			blocks.append([ founders[i], int(new_starts[i]), int(new_ends[i]) ])

		return chrom, blocks

	# core function of this module: parse the hap file passed to the class constructor
	def parse(self):
		# initialise empty dictionary of haplotype blocks:
		blocks = { 0: Genome(self._sample), 1: Genome(self._sample) }
		# keep list of chromosomes visited already, to detect phasing
		chroms_seen = []

		# loop on hapfile and parse informative lines (skipping display-only lines)
		for line in self.file:
			if not self._is_chromline(line):
				continue
			else:
				chrom, newblocks = self._parse_chromline(line)
				phase = self._get_phase( chrom, chroms_seen )
				chroms_seen.append( chrom )
				for b in newblocks:
					blocks[phase].add_segment(chrom, b[1], b[2], phase, b[0])

		self.hapblocks = blocks

	def get_recombs(self):
		for phase in self.hapblocks:
			for recomb in self.hapblocks[phase].get_recombinations():
				yield recomb

	def get_hapblocks(self):
		for phase in self.hapblocks:
			for block in self.hapblocks[phase]:
				yield block
