
import TransposableElements as TE
import argparse
import gzip
import pyranges as pr
import sys
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("rmsk_file", type=str, help="RepeatMasker reference file")
parser.add_argument("meta_file", type=str, help="Tab-separated file with two columns: element names, and their 'meta' types (e.g., 'internal' and 'LTR')")
parser.add_argument("-g", "--max_gap", type=int, help="Maximum gap [bp] between elements to be merged")
args = parser.parse_args()

if args.rmsk_file.endswith(".gz"):
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
	lineno = 0
	for line in inf.readlines():
		lineno += 1
		#if lineno % 2000 == 0:
		#	print("Read {} subelements".format(lineno))
		entry = line.rstrip('\n').split('\t')
		elem_ID = int(entry[9]) 
		
		try:
			entry_elem = TE.ERV(TE.Subelement(TE.ElementType(entry[5], entry[6], entry[4], meta_types[entry[4]]), TE.GenomicPosition(entry[0], int(entry[1]), int(entry[2]), entry[3]), [int(entry[7]), int(entry[8])]), elem_ID)
			
			if elem_ID in elements:
				elements[elem_ID].merge(entry_elem)
			else:
				elements[elem_ID] = entry_elem
		except Exception as e:
			print("Error on line {}".format(lineno), file=sys.stderr)
			print(e, file=sys.stderr)
			sys.exit(1)

ERV_elem_gp = [erv.span() for erv in elements.values()]

erv_df = pd.DataFrame({
	'Chromosome': [e.chrom for e in ERV_elem_gp],
	'Start': [e.start for e in ERV_elem_gp],
	'End': [e.end for e in ERV_elem_gp],
	'Strand': [e.strand for e in ERV_elem_gp],
	'ID': [e.id for e in elements.values()],
	'Struct': [e.meta_str() for e in elements.values()]
})
ERVs = pr.PyRanges(erv_df)

def update_pr(changed_id, removed_id):
	global ERVs
	#print("{}\t{}".format(changed_id, removed_id))
	new_elem = elements[changed_id].span().pr()
	new_elem.ID = changed_id
	new_elem.Struct = elements[changed_id].meta_str()
	#print(new_elem)
	print("Merging {} into {}".format(removed_id, changed_id))
	ERVs = pr.concat([pr.PyRanges(ERVs.df.loc[~ERVs.df['ID'].isin([changed_id, removed_id])]), new_elem])

incomplete_IDs = [e.id for e in elements.values() if not e.is_complete()]
print("Incomplete elems: {}".format(len(incomplete_IDs)))
processed = 0
for elem_ID in incomplete_IDs:
	processed += 1
	if processed % 1000 == 0:
		print("Processing {}".format(processed))
	if elem_ID not in elements:
		continue
	e = elements[elem_ID]
	
	if e.missing_5prime_LTR():
		adjacent_elem_IDs = ERVs.overlap(e.span().flanked_5(args.max_gap).pr(), strandedness="same").ID
		if e.strand == "-":
			adjacent_elem_IDs = adjacent_elem_IDs[::-1]
		adjacent_elems = [elements[i] for i in adjacent_elem_IDs if i != e.id and e.compatible_with(elements[i])]
		complementary_elems = [a for a in adjacent_elems if a.missing_3prime_LTR() and not a.missing_5prime_LTR()]
		if len(complementary_elems) > 0:
			e.merge(complementary_elems[0])
			del elements[complementary_elems[0].id]
			update_pr(elem_ID, complementary_elems[0].id)
		else:
			solo_elems = [a for a in adjacent_elems if a.is_solo_LTR()]
			if len(solo_elems) > 0:
				e.merge(solo_elems[0])
				del elements[solo_elems[0].id]
				update_pr(elem_ID, solo_elems[0].id)
			else:
				full_elems = [a for a in adjacent_elems if a.is_fully_structured()]
				if len(full_elems) > 0:
					e.merge(full_elems[0])
					del elements[full_elems[0].id]
					update_pr(elem_ID, full_elems[0].id)
					
	if e.missing_3prime_LTR():
		adjacent_elem_IDs = ERVs.overlap(e.span().flanked_3(args.max_gap).pr(), strandedness="same").ID
		if e.strand == "-":
			adjacent_elem_IDs = adjacent_elem_IDs[::-1]
		adjacent_elems = [elements[i] for i in adjacent_elem_IDs if i != e.id and e.compatible_with(elements[i])]
		complementary_elems = [a for a in adjacent_elems if a.missing_5prime_LTR() and not a.missing_3prime_LTR()]
		if len(complementary_elems) > 0:
			e.merge(complementary_elems[0])
			del elements[complementary_elems[0].id]
			update_pr(elem_ID, complementary_elems[0].id)
		else:
			solo_elems = [a for a in adjacent_elems if a.is_solo_LTR()]
			if len(solo_elems) > 0:
				e.merge(solo_elems[0])
				del elements[solo_elems[0].id]
				update_pr(elem_ID, solo_elems[0].id)
			else:
				full_elems = [a for a in adjacent_elems if a.is_fully_structured()]
				if len(full_elems) > 0:
					e.merge(full_elems[0])
					del elements[full_elems[0].id]
					update_pr(elem_ID, full_elems[0].id)

for elem in elements.values():
	subs = elem.sub
	if elem.strand == "-":
		subs = elem.sub[::-1]
	for subelem in subs:
		print("\t".join([str(x) for x in subelem.pos, sub.type.Name, sub.type.Class, subelem.type.Family, subelem.rep_start, subelem.rep_end, elem.id]))
