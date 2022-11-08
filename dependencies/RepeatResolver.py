## the purpose of this script is to take as an input the .tsv file produced by USEARCH_local and parse some output .fastas

## after USEARCH_local has purformed local alignment of a MARS circularly permuted cluster centroid (query) against a larger concatemer (target)
## where the centroid can be defined as "unit length" (UL) and the concatemer as "over unit length" (>UL)
## the output .tsv is defined by these settings: '-userfields target+id+tlo+thi+trow+ql+tl'
	## so per line (hit), we have a tab seperated list of:
		## the seqID of the >UL, %ID, start coord on hit, end coord on hit, hit sequence + INDELS, query starting length, target starting length
	## so the >UL length and adjusted coverage can be extracted from the first field

## note, because of how the MARS+USEACH_local is being run, there will always only be one target sequence and hence one set of hits

## the script will mostly work off the start/end coordinates (compared to the >UL length) to try and determine if the alignment makes sense
	## the point being that if these values dont make sense, that the >UL was actually somehow wrong and it should be omitted from the data
	## otherwise, the >UL will be sub-divided into ~full length and ~fragment length .fastas for individual output

## flow:
	## 0a) check that the hits from USEARCH_local are continuous and not overlapping beyond some small leeway (iff = proceed)
		## the hits can't be too overlapped or too strongly spaced
	## 0b) count how many unit-length hits are expected given the query and target lengths + X% error
		## check that the number of unit-length hits (+error) satisfy this expectation (iff = proceed)
	## 1) scrubbing the hits for INDELs ('-'), save the unit lengths to a .fasta file
		## being sure to rename the seqID to account for the change in length and hence AC (+ a new identifier?)
	## 2) also save the remaining sub-unit-lengths too to their own .fastas

## note, for all output .fastas, the adjusted coverage will be computed from:
	## the product of the intial adjusted coverage and the ratio of the initial target length to the mean unit-length length
	## so fragments and unit length hits will have the same adjusted coverage

import argparse
import math
from statistics import mean

def main(input_tsv_name, output_fasta_name_prefix, CCC_pair_identifier, overlap_spacing_error_fxn, full_len_error_fxn, debug):
	
	## declare variables
	
	##input_tsv_name = 'CCCpair_sense_0_USL.tsv'
	##output_fasta_name_prefix = 'SRR11060619_sub_resolved_cluster_0'
	##CCC_pair_identifier = '0'
	##overlap_spacing_error_fxn = 0.05
	##full_len_error_fxn = 0.1
	##debug = 1
	
	overlap_percentage = round((overlap_spacing_error_fxn*100),0)
	
	output_fasta_name_prefix = output_fasta_name_prefix + "_pair_" + CCC_pair_identifier + "_subID_"
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry.split('\t')
		##
		return open_list
	##
	
	input_tsv_list = opener(input_tsv_name)
	
	## determine the >UL length and starting coverage, and query length
	
	target_start_len = int(input_tsv_list[0][6])
	target_start_cov = float(input_tsv_list[0][0].split('_')[7])
	
	query_len = int(input_tsv_list[0][5])
	
	overlap_nt_len = int(query_len*overlap_spacing_error_fxn)
	
	## first sort the start-end co-ords in order (keeping the pairing)

	coords_list = []
	temp_list = []
	for hit in input_tsv_list:
		start_coord = int(hit[2])
		end_coord = int(hit[3])
		temp_list.append(start_coord)
		temp_list.append(end_coord)
		coords_list.append(temp_list)
		temp_list = []
	##
	
	coords_list = sorted(coords_list, key=lambda x: x[0])
	
	def overlapcheck(tsv_list_in, coords_list_in, check_nt_len):
		passing_overlap = True
		for coord_pair_ind, coord_pair in enumerate(coords_list_in):
			current_hit_end = coord_pair[1]
			if (coord_pair_ind == len(coords_list_in)-1):
				next_hit_start = target_start_len
			else:
				next_hit_start = coords_list_in[coord_pair_ind+1][0]
			##
			overlap = abs(next_hit_start-current_hit_end)
			if (overlap > check_nt_len):
				passing_overlap = False
				if debug:
					print("hits did not satisfy the overlap/spacing criteria")
					print("found overlap/spacing of: " + str(overlap) + "nt")
				##
				break
		##
		if passing_overlap:
			if debug:
				print("hits did satisfy the overlap/spacing criteria")
		##
		return passing_overlap
	##
	
	def ULcountcheck(coords_list_in, target_len_in, query_len_in, error_fxn):
		passing_UL_count = True
		min_query_len = round(query_len_in*(1-error_fxn))
		max_query_len = round(query_len_in*(1+error_fxn))
		min_expected_UL_count = math.floor(target_len_in/max_query_len)
		max_expected_UL_count = math.floor(target_len_in/min_query_len)
		hit_len_list = []
		if debug:
			print("expecting at least " + str(min_expected_UL_count) + " and at most " + str(max_expected_UL_count) + " unit-length hits")
		##
		unit_length_hit_count = 0
		for coord_pair in coords_list:
			hit_len = coord_pair[1] - coord_pair[0]
			if (hit_len >= min_query_len) and (hit_len <= max_query_len):
				unit_length_hit_count += 1
				hit_len_list.append(hit_len)
		##
		if not (unit_length_hit_count >= min_expected_UL_count) and (unit_length_hit_count <= max_expected_UL_count):
			passing_UL_count = False
			if debug:
				print("did not find the expected number of unit-length hits: " + str(unit_length_hit_count))
		##
		mean_len = 1
		if hit_len_list:
			mean_len = round(mean(hit_len_list))
		##
		return passing_UL_count, mean_len
	##

	## determine if the hits are sufficiently non-overlapping and continuous
	## and then check if the hits contain the expected number of unit length (~query) hits
	
	if debug:
		print("-----=====-----")
		print("determining if the hits are continuous")
		print("allowing for overlap/spacing within +/- " + str(overlap_percentage) + "% unit length (" + str(overlap_nt_len) + "nt)")
	##

	if overlapcheck(input_tsv_list, coords_list, overlap_nt_len):
		if debug:
			print("-----=====-----")
			print("determining if the hits contain the expected number of ~unit-length hits")
		##
		hits_passed, mean_UL_len = ULcountcheck(coords_list, target_start_len, query_len, full_len_error_fxn)
		##
		if hits_passed:
			if debug:
				print("found the expected number of unit-length hits")
				print("-----=====-----")
				print("extracting hits into unit-length and fragment-length .fastas")
		else:
			if debug:
				print("did not the expected number of unit-length hits")
				print("no new .fastas will be saved")
				mean_UL_len = 1
				input_tsv_list = []
		##
	##
	else:
		print("no new .fastas will be saved")
		mean_UL_len = 1
		input_tsv_list = []
	##
	
	## now let's extract the hits
	## let's assume that if hits were found, that the adjusted coverage should be computed from the mean hit length for all the output sequences
	
	new_AC = round(target_start_cov * (target_start_len/mean_UL_len),4)
	
	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##
	
	def hitextractor(tsv_list_in, query_len_in, error_fxn, input_AC, input_prefix):
		output_fasta_list = []
		seq_ID_counter = 0
		digit_count = len(str(len(tsv_list_in)))
		for hit in tsv_list_in:
			seq_ID_str = str(seq_ID_counter).zfill(digit_count)
			init_seqID_start = '_'.join(hit[0].split('_')[:5]) + "_" ## extract the front of the seqID until the old length
			init_seqID_end = "_" + '_'.join(hit[0].split('_')[-4:]) + "_subID_" ## extract the end of the seqID after the old coverage
			hit_seq = hit[4].replace('-', '')
			len_hit_seq = len(hit_seq)
			if (len_hit_seq < round(query_len_in*(1-error_fxn))):
				suffix = "_FR"
			else:
				suffix = "_UL"
			##
			suffix = seq_ID_str + suffix
			new_seqID = ">" + init_seqID_start + str(len_hit_seq) + "_AC_" + str(input_AC) + init_seqID_end + suffix
			output_fasta_list.append(new_seqID)
			output_fasta_list.append(hit_seq)
			output_fasta_filename = input_prefix + suffix + ".fasta"
			saver(output_fasta_filename, output_fasta_list)
			output_fasta_list = []
			seq_ID_counter += 1
		##
		print("saved " + str(seq_ID_counter) + " .fasta files")
	##
	
	hitextractor(input_tsv_list, query_len, full_len_error_fxn, new_AC, output_fasta_name_prefix)
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-input_tsv', type=str, help='input USEARCH_local derived .tsv e.g.: "CCCpair_sense_0_USL.tsv"')
	parser.add_argument('-prefix', type=str, default = 'output_seq', help='prefix string for all output .fasta filenames - default = "output_seq"')
	parser.add_argument('-pair_ID', type=str, default = 'X', help='string identifier for the specific C-CC pair assessed - default = "X"')
	parser.add_argument('-overlap_fxn', type=float, default = 0.05, help='fraction of unit-length length permitted for hit overlap/spacing - default = 0.05 = 5 percent')
	parser.add_argument('-unit_length_fxn', type=float, default = 0.1, help='fractional boundaries around the unit-length length permitted for unit-length hit identification - default = 0.1 = 10 percent')
	parser.add_argument('-debug', type=int, default = 1, help='print more detailed results - default = 1')
	
	args = parser.parse_args()
	main(args.input_tsv, args.prefix, args.pair_ID, args.overlap_fxn, args.unit_length_fxn, args.debug)
##















































