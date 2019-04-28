
import TransposableElements as TE
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("rmsk_file", type="str", help="RepeatMasker reference file")
parser.add_argument("meta_file", type-"str", help="Tab-separated file with two columns: element names, and their 'meta' types (e.g., 'internal' and 'LTR'")
args = parser.parse_args()

if rmsk_file.endswith(".gz"):
	fopen = gzip.open
else:
	fopen = open

meta_types = {}
with open(args.meta_file) as mf:
	for line in mf:
		(name, meta) = line.strip('\n').split('\t')
		meta_types[name] = meta

# store all elements in a dict by chrom and element ID
elements = {} 

# read through reference file and store TEs 
with fopen(args.rmsk_file, 'rt') as inf:
	for line in inf.readlines():
		entry = line.rstrip('\n').split('\t')

		# parse entry

		entry_elem = TE.TransposableElement(Subelement(), elem_ID)
		
		if not chrom in elements:
			elements[chrom] = {}
		
		if elem_ID in elements[chrom]:
			elements[chrom][elem_ID] = elements[chrom][elem_ID].merge(entry_elem)
		else:
			elements[chrom][elem_ID] = entry_elem

# iterate through incomplete elements
for chrom in elements.keys():
	pass
