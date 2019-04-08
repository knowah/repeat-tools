#!/usr/bin/python3

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="Input (tab-separated) RepeatMasker file (must be sorted)")
parser.add_argument("breaks",  type=int, default=1000000, help="Breakpoints to merge upon")
parser.add_argument("-s", "--seq_gap", type=int, default=0, help="Maximum (reference sequence) gap in bp between broken subelements")
parser.add_argument("-r", "--rep_gap", type=int, default=1, help="Maximum gap between repEnd & repStart of broken subelements")
parser.add_argument("-f", "--fix", action="store_true", help="Merge subelements at the breakpoint when unambiguous [NOT IMPLEMENTED]")
args = parser.parse_args()

if args.fix:
	raise InputError("--fix not implemented! re-run without this option")

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

with open(args.infile) as inf:
	prev = process_line(inf.readline())
	merge_elements = False
	elem_ids = (None, None)
	for line in inf:
		curr = process_line(line)
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

		if merge_elements and curr.element_ID == elem_ids[1]:
			curr.element_ID = elem_ids[0]
		else:
			merge_elements = False
		
		print(prev)
		prev = curr
	
	print(prev) # print last line
