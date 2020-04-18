#!/usr/bin/env python3

import argparse
import sys
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("infile",  type=str, help="Input (tab-separated) RepeatMasker file (must be sorted)")
parser.add_argument("breaks",  type=int, default=500000, help="Breakpoints to merge upon [default: 500000]")
parser.add_argument("-s", "--seq_gap", type=int, default=0, help="Maximum (reference sequence) gap in bp between broken subelements")
parser.add_argument("-r", "--rep_gap", type=int, default=1, help="Maximum gap between repEnd & repStart of broken subelements")
parser.add_argument("-f", "--fix", action="store_true", help="Merge subelements at the breakpoint when unambiguous [NOT IMPLEMENTED]")
parser.add_argument("--alt_format", action="store_true", help="Use alternative input format (e.g. from mouse strain genomes)")
args = parser.parse_args()

if args.fix:
	raise NotImplementedError("--fix is not yet implemented! Re-run without this option.")

FIELDS = ['chrom', 'start', 'end', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'element_ID']
class RM_Entry(object):
	__slots__ = FIELDS
	def __str__(self):
		return "\t".join([str(self.__getattribute__(x)) for x in FIELDS])

def process_line(ln):
	tokens = ln.rstrip('\n').split('\t')
	# chrom start end strand repName repClass repFamily repStart repEnd element.ID
	entry = RM_Entry()
	entry.chrom  = tokens[5]
	entry.start  = int(tokens[6])
	entry.end    = int(tokens[7])
	entry.strand = tokens[9]
	entry.repName   = tokens[10]
	entry.repClass  = tokens[11]
	entry.repFamily = tokens[12]
	entry.repStart = int(tokens[13])
	entry.repEnd   = int(tokens[14])
	entry.element_ID = int(tokens[16])
	return entry
	#return [(tokens[i] if i not in num_fields else int(tokens[i])) for i in [5,6,7,9,10,11,12,13,14,16]]

def process_line_alt(ln):
	tokens = ln.rstrip('\n').split('\t')
	#chrom chromStart chromEnd name score strand swScore milliDiv milliDel milliIns genoLeft repClass repFamily repStart repEnd repLeft
	entry = RM_Entry()
	entry.chrom  = tokens[0]
	entry.start  = int(tokens[1])
	entry.end    = int(tokens[2])
	entry.strand = tokens[5]
	entry.repName   = tokens[3]
	entry.repClass  = tokens[11]
	entry.repFamily = tokens[12]
	entry.repStart = abs(int(tokens[13]))
	entry.repEnd   = int(tokens[14])
	entry.element_ID = lineno
	return entry
		

# process RepeatMasker file by line
if args.infile.endswith(".gz"):
	fopen = gzip.open
else:
	fopen = open

read_entry = process_line if not args.alt_format else process_line_alt

with fopen(args.infile, 'rt') as inf:
	lineno = 1
	firstln = inf.readline()
	if firstln.startswith('#'):
		firstln = inf.readline()
	prev = read_entry(firstln)
	merge_elements = False
	elem_ids = (None, None) # when merge_elements == True, this stores the element IDs to merge

	for line in inf:
		lineno += 1
		curr = read_entry(line)

		# fail if input file not sorted
		if prev.chrom == curr.chrom and curr.start < prev.start:
			raise IOError("ERROR: Entries out of order (line {})".format(lineno))

		# check if this element and the previous one need to be merged
		if not merge_elements \
		   and curr.start % args.breaks <= args.seq_gap \
		   and curr.chrom == prev.chrom \
		   and curr.strand == prev.strand \
		   and curr.start - prev.end <= args.seq_gap \
		   and curr.repName == prev.repName \
		   and curr.repClass == prev.repClass \
		   and curr.repFamily == prev.repFamily \
		   and ((curr.strand == "+" and curr.repStart - prev.repEnd <= args.rep_gap) or \
		        (curr.strand == "-" and prev.repStart - curr.repEnd <= args.rep_gap)):
			merge_elements = True
			elem_ids = (prev.element_ID, curr.element_ID)
		
		# merge elements if the element ID matches the ID to be changed
		if merge_elements and curr.element_ID == elem_ids[1]:
			curr.element_ID = elem_ids[0]
		else:
			merge_elements = False
		
		print(prev) # output new entry for the previous line
		prev = curr
	
	print(prev) # print last line
