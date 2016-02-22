#! /usr/bin/env python
"""
genome.py
Classes for representing genomes from multiparent populations as collections of founder-labelled segments.
"""

from __future__ import print_function

import os
import sys
import warnings
from collections import defaultdict

from bx.intervals import *
#from intervaltree import *

SEX_UNKNOWN = 0
SEX_MALE = 1
SEX_FEMALE = 2

PHASE_UNKNOWN = -1
PHASE_MATERNAL = 0
PHASE_PATERNAL = 1

AUTOSOMES = [ "chr"+str(c) for c in range(1,20) ]
CHROMS = AUTOSOMES + [ "chrX","chrY","chrM" ]
PLOIDY = {	SEX_UNKNOWN:	{ "chrX": 2, "chrY": 1, "chrM": 1 },
			SEX_MALE:		{ "chrX": 1, "chrY": 0, "chrM": 1 },
			SEX_FEMALE:		{ "chrX": 1, "chrY": 0, "chrM": 1 } }

FOUNDERS = { 	"A": "A/J", "B": "C57BL/6J", "C": "129S1/SvImJ",
				"D": "NOD/ShiLtJ", "E": "NZO/HlLtJ",
				"F": "CAST/EiJ", "G": "PWK/PhJ", "H": "WSB/EiJ" }

class HaplotypeSegment(Interval):

	def __init__(self, start, end, **kwds):
		super(HaplotypeSegment, self).__init__(start, end, **kwds)
		self.phase = None
		self.chrom = None
		self.founder = None
		self.attributes = {}

	def __repr__(self):
		rez = "HaplotypeSegment( {}:{}-{} [{}:{}] )".format(self.chrom, self.start, self.end, self.phase, self.founder)
		return rez

	def __str__(self):
		return self.__repr__()

class Genome:

	def __init__(self, iid, sex = SEX_UNKNOWN):
		self.iid = iid
		self.sex = sex
		self._blocks = { c: IntervalTree() for c in CHROMS }

	def __repr__(self):
		rez = "Genome of individual '{}'\n\n".format(self.iid)
		rez += "Haplotype segments\n---------------------\n"
		for chrom in CHROMS:
			if chrom in self._blocks:
				rez += "  {}:".format(chrom)
				nb = 0
				for b in self._blocks[chrom].find(0,1e9):
					rez += "\t      {}:({}, {})\n".format(b.founder, b.start, b.end)
					nb += 1
				if not nb:
					rez += "\n"
		rez += "\n\n"
		return rez

	def __iter__(self):
		for chrom in CHROMS:
			if chrom in self._blocks:
				for b in self._blocks[chrom].find(0, 1e9):
					yield b

	def _get_founder_name(self, query):
		if query in FOUNDERS:
			return query
		elif query in FOUNDERS.values():
			founder = None
			choices = FOUNDERS.keys()
			i = 0
			while not founder:
				if query == FOUNDERS[ choices[i] ]:
					founder = choices[i]
				i += 1
			return founder
		else:
			return None

	def _check_segment(self, chrom, start, end):
		query_ok = True
		query_ok &= chrom in CHROMS
		query_ok &= start >= 0
		query_ok &= end >= start
		return query_ok

	def add_segment(self, chrom, start = None, end = None, phase = 0, founder = None, extras = {}, check = True):
		if not self._check_segment(chrom, start, end):
			raise ValueError("Nonsensical interval.")
		if (not founder in FOUNDERS) and (founder is not None) and check:
			raise ValueError("Unknown founder: {}.".format(founder))
		#new_segment = HaplotypeSegment(start, end, chrom = chrom, phase = phase, founder = founder)
		new_segment = HaplotypeSegment(start, end)
		new_segment.chrom = chrom
		new_segment.founder = founder
		new_segment.phase = phase
		print(str(new_segment))
		self._blocks[chrom].add_interval(new_segment)

	def insert_segment(new_segment):
		self._blocks[new_segment.chrom].add_interval(new_segment)

	def get_segments(self, chrom, start = 0, end = 1e9):
		if isinstance(chrom, Interval):
			start = chrom.start
			end = chrom.end
			chrom = chrom.chrom
		if end < start:
			x = start
			start = end
			end = x
		if self._check_segment(chrom, start, end):
			olaps = [ o for o in self._blocks[chrom].find(start, end) ]
			expect_ploidy = 2
			if chrom not in AUTOSOMES:
				expect_ploidy = PLOIDY[self.sex][chrom]
			#if expect_ploidy != len(olaps):
			#	warnings.warn("Number of overlapping intervals ({}) is less than expected for this chromosome {} and sex {}.".format(len(olaps), chrom, self.sex), RuntimeWarning)
			return olaps
		else:
			raise ValueError("Nonsensical query: check chromosome name and verify that start > 0 and end >= start.")

	def get_recombinations(self):
		chrom = [ c for c in CHROMS if c in self._blocks ]
		for c in chrom:
			last = [ None, None ]
			for block in self.get_segments(c):
				if block.phase is None:
					continue
				if last[block.phase]:
					recomb = HaplotypeSegment(last[block.phase].end, block.start)
					recomb.chrom = c
					recomb.phase = block.phase
					recomb.attributes = { "from": last[block.phase].founder, "to": block.founder }
					yield recomb
				last[block.phase] = block

	def get_genotype_at(self, chrom, start = None):
		blocks = self.get_segments(chrom, start, start)
		ploidy = 2
		if chrom not in AUTOSOMES:
			ploidy = PLOIDY[self.sex][chrom]
		gty = [None]*ploidy
		curr_idx = 0
		if ploidy == 2:
			for b in blocks:
				if b.phase == 0:
					gty[0] = "{}".format(b.founder)
				elif b.phase == 1:
					gty[1] = "{}".format(b.founder)
				elif curr_idx < len(gty):
					gty[curr_idx] = "{}".format(b.founder)
					curr_idx += 1
		else:
			for i,b in enumerate(blocks):
				gty[i] = "{}".format(b.founder)
		return gty


if __name__ == "__main__":

	mouse = Genome("CC001/Unc", sex = SEX_UNKNOWN)
	mouse.add_segment("chr1",0,100000,0,"A")
	mouse.add_segment("chr1",100000,200000,0,"C")
	mouse.add_segment("chr1",0,100000,1,"B")
	print(str(mouse))

	for block in mouse:
		print(block)

	print(mouse.get_genotype_at("chr1", 500))
