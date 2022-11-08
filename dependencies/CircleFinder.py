## the purpose of this script is to make use of the 3' k-mer repeat "artifact" that rnaSPAdes generates on circular/tandem-repeat contigs
	## the last k-len that spades.log outputs is the k-len used to assemble the contigs
	## if the de Bruijn graph was circular, than the resulting contig has the first k nucleotides repeated at the 3' end perfectly
		## so, by looking for an end repeat of the terminal 3' k-mer at the 5' end, one can enrich for circular contigs

## "simple circular contig":

##		<-------------unit_length------------->
##		START----------------------------------|TERM

## in this case, for each contig, take the k-mer "TERM" (where k-len is from spades.log) and see if it matches "START"
	## then trim off TERM, yielding a unit-length contig
	## so in this example TERM==START


## "extended simple circular contig":

##		<-------------unit_length------------->
##		EXT-START-----------------------------|EXT-TERM

## in this case, for whatever reason, a quasi-circle graph was made (?) and TERM is embedded within the 5' end vs at the 5' end
	## here, the task would be to scan TERM along the contig until a match is found (START) - that is not the 3' TERM
	## then to look 5' of the match to identify the "extention" (EXT), append that to TERM, and trim EXT-TERM off the 3' end
		## yeilding a unit-length contig
	## so then the "simple" case is the scenario where len(EXT)=0


## "tandem circular contig"

##		<-------------unit_length-------------><-------------unit_length------------->
##		EXT-START-----------------------------|EXT-START-----------------------------|EXT-TERM
##		...............repeat_1...............|...............repeat_2...............|

## in this case, for whatever reason, say erroneous concatenation of quasi-species, the contig is multiple unit_lengths long
	## so repeat_1 and repeat_2 contain some polymorphisms that prevented rnaSPAdes from collapsing them
	## in this case once EXT-START is found, the whole contig should be scanned for instances of EXT-START
	## then, for each instance of EXT-START, a unit_length monomer can be extracted by iteratively making note of the coordinates of each monomer
		## the expectation then would be that each monomer is roughly the same length - this would be as a check

## lastly, given that the concatenation of quasispecies relies on polymorphisms
## it stands to reason that there could be polymorphisms between (EXT-)TERM and (EXT-)START
	## so, it might make sense to search for these k-mers using a fuzzy string matching approach - say fuzzywuzzy or the python regex package

## so general approach:
	## 0) load the rnaSPAdes contigs:
		## check (and correct?) for singleline multifasta format
	## 1) set variables:
		## k-len
		## to only search for "simple" circles
		## or to allow for "extensions"
		## to then filter these hits and search for search for "tandem" circles
			## when "tandem" is chosen, "simple" circles will also be found
			## enforce a length check of the monomers from "tandem" circles
				## allow the fractional length check threshold to be set (default 5%?)
		## to allow for fuzzy k-mer matching (+ what degree of fuzziness) - this is experimental and put in a seperate version of this script
	## 2) given the settings, search for the circular contigs
	## 3) rename the identified circles to contain their inferred unit length (UL)
		## and adjust the coverage (AC) value by dividing by the ratio of UL to original length
		## append specified name strings to the new seqIDs
	## 4) save the indentified circles into one final .fasta file and report the stats of how many of each kind of circle were found

## CircleFinder_2 will have the added functionality of being able produce a file of the determind "linear" contigs left in the input
	## this will be an optional setting
	## these contigs will undergo similar renaming
	## this will be done by appending to the simplecircle and extendedcircle functions a recorder of the indices that were determined "circular"
	## then a dedicated function will keep the indices that weren't circular and rename them

## import libraries:

from statistics import mean
import copy
import argparse
	
def main(input_fasta_name, k_len, simple_circ, tandem_circ, tandem_thr, min_len, keep_linear, new_seqID_prefix, output_circ_fasta_name, output_lin_fasta_name, debug, verb_debug):
	
	## define variables:
	
	##input_fasta_name = 'SRR11060619_sub_rnaSPAdes.fasta'
	##k_len = 73
	##simple_circ = 0
	##tandem_circ = 1
	##tandem_thr = 0.1
	##min_len = 100
	##keep_linear = 1
	##new_seqID_prefix = 'SRR11060619_sub'
	##output_circ_fasta_name = 'SRR11060619_sub_rnaSPAdes_circles.fasta'
	##output_lin_fasta_name = 'SRR11060619_sub_rnaSPAdes_linears.fasta'
	##debug = 1
	##verb_debug = 0
	
	## open file
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry
		##
		return open_list
	##
	
	input_fasta_list = opener(input_fasta_name)
	
	## convert to singleline (based on https://stackoverflow.com/a/50856787):
	
	def singleline(in_fasta_list):
		out_fasta_list = []
		new_seq_line = []
		for line in in_fasta_list:
			if line.startswith('>'):
				if new_seq_line:
					out_fasta_list.append(''.join(new_seq_line))
					new_seq_line = []
				out_fasta_list.append(line)
			else:
				new_seq_line.append(line.strip())
			##
		if new_seq_line:
			out_fasta_list.append(''.join(new_seq_line))
		##
		return out_fasta_list
	##
	
	input_fasta_list = singleline(input_fasta_list)
	
	## function for searching for "simple" circular contigs:
	
	def simplecircle(in_fasta_list):
		out_simcirc_list = []
		out_simcirc_ind_list = []
		for line_ind, line in enumerate(in_fasta_list):
			if (line_ind%2):
				end_kmer = line[-k_len:]
				start_kmer = line[:k_len]
				if (end_kmer == start_kmer):
					trimmed_seq = line[:-k_len]
					out_simcirc_list.append(in_fasta_list[line_ind-1])
					out_simcirc_list.append(trimmed_seq)
					out_simcirc_list.append(end_kmer)
					out_simcirc_ind_list.append(line_ind-1)
					out_simcirc_ind_list.append(line_ind)
					if verb_debug:
						print("-----=====-----")
						print("simple circle found")
						print("end k-mer   = " + end_kmer)
						print("start k-mer = " + start_kmer)
						print("seqID = " + str(in_fasta_list[line_ind-1]))
						print("input seq  = " + str(line))
						print("output seq = " + str(trimmed_seq))
		##
		return out_simcirc_list, out_simcirc_ind_list
	##
	
	## function for searching for "extended" circular contigs:
	
	def extendedcircle(in_fasta_list):
		out_extcirc_list = []
		out_extcirc_ind_list = []
		first_ext_found = False
		for line_ind, line in enumerate(input_fasta_list):
			if (line_ind%2):
				init_end_kmer = line[-k_len:]
				for scan_kmer_ind in range(len(line)-k_len+1):
					scan_kmer = line[scan_kmer_ind:scan_kmer_ind+k_len]
					if (scan_kmer == init_end_kmer):
						if (scan_kmer_ind != (len(line)-k_len)):
							if not first_ext_found:
								first_ext_found = True
								ext_kmer = line[:scan_kmer_ind]
								ext_end_kmer = ext_kmer + init_end_kmer
								new_end_kmer = line[-len(ext_end_kmer):]
								if (ext_end_kmer == new_end_kmer):
									trimmed_seq = line[:-len(new_end_kmer)]
									out_extcirc_list.append(in_fasta_list[line_ind-1])
									out_extcirc_list.append(trimmed_seq)
									out_extcirc_list.append(new_end_kmer)
									out_extcirc_ind_list.append(line_ind-1)
									out_extcirc_ind_list.append(line_ind)
									if verb_debug:
										print("-----=====-----")
										print("extended circle found")
										if not scan_kmer_ind:
											print("simple (non-extended) case")
										print("seqID = " + str(input_fasta_list[line_ind-1]))
										print("initial end k-mer = " + init_end_kmer)
										print("matchd scan k-mer = " + scan_kmer)
										print("input seq   = " + str(line))
										print("matchd scan k-mer start index = " + str(scan_kmer_ind))
										print("extended kmer seq = " + ext_kmer)
										print("exte-end kmer seq = " + ext_end_kmer)
										print("new end kmer seq  = " + new_end_kmer)
										print("trimmed seq = " + trimmed_seq)
				first_ext_found = False
		##
		return out_extcirc_list, out_extcirc_ind_list
	##
	
	numb_input_seq = str(int(len(input_fasta_list)/2))
	if simple_circ:
		circ_list, circ_ind_list = simplecircle(input_fasta_list)
		if debug:
			print("-----=====-----")
			print("running simple 'circle' discovery on " + numb_input_seq + " input sequences")
			numb_output_seq = str(int(len(circ_list)/3))
			print("found " + numb_output_seq + " simple 'circles'")
	else:
		circ_list, circ_ind_list = extendedcircle(input_fasta_list)
		if debug:
			print("-----=====-----")
			print("running extended 'circle' discovery on " + numb_input_seq + " input sequences")
			numb_output_seq = str(int(len(circ_list)/3))
			print("found " + numb_output_seq + " simple or extended 'circles'")
	##
	
	if not circ_list:
		print("no 'cirlces' identified")
		print("qutting")
		quit()
	##

	## function for searching for "tandem" circular contigs from the identified 'circles':
		## input is a list of original_seqID, trimmed_seq, and end_repeat_kmer
	
	def tandemcircle(in_circ_list):
		out_no_tandcirc_list = []
		out_yes_tandcirc_list = []
		for line_ind, line in enumerate(circ_list):
			if line.startswith('>'):
				seq_line = circ_list[line_ind+1]
				end_repeat_kmer = circ_list[line_ind+2]
				ER_kmer_len = len(end_repeat_kmer)
				matched_scan_kmer_count = 0
				matched_scan_kmer_inds_list = []
				for scan_kmer_ind in range(len(seq_line)-ER_kmer_len+1):
					scan_kmer = seq_line[scan_kmer_ind:scan_kmer_ind+ER_kmer_len]
					if (scan_kmer == end_repeat_kmer):
						matched_scan_kmer_count += 1
						matched_scan_kmer_inds_list.append(scan_kmer_ind)
				if (matched_scan_kmer_count > 1):
					temp_yes_tandcirc_list = []
					temp_yes_tandcirc_len_list = []
					disimilar_len = False
					matched_scan_kmer_inds_list.append(len(seq_line))
					for junct_kmer_pos_ind, junct_kmer_pos in enumerate(matched_scan_kmer_inds_list[:-1]):
						start_ind = junct_kmer_pos
						end_ind = matched_scan_kmer_inds_list[junct_kmer_pos_ind+1]
						monomer_seq = seq_line[start_ind:end_ind]
						temp_yes_tandcirc_list.append(line)
						temp_yes_tandcirc_list.append(monomer_seq)
						temp_yes_tandcirc_list.append(end_repeat_kmer)
						temp_yes_tandcirc_len_list.append(len(monomer_seq))
					norm_monomer_len_list = [x / temp_yes_tandcirc_len_list[0] for x in temp_yes_tandcirc_len_list]
					for norm_len in norm_monomer_len_list:
						if (norm_len <= (1-tandem_thr)) or (norm_len >= (1+tandem_thr)):
							disimilar_len = True
					if not disimilar_len:
						out_yes_tandcirc_list.extend(temp_yes_tandcirc_list)
						if verb_debug:
							monomer_count = len(temp_yes_tandcirc_len_list)
							monomer_mean_len = round(mean(temp_yes_tandcirc_len_list),0)
							print("-----=====-----")
							print("tandem circle found")
							print("seqID = " + str(line))
							print("input seq        = " + str(seq_line))
							print("junction kmer    = " + end_repeat_kmer)
							print("last monomer seq = " + str(monomer_seq))
							print("monomer mean length = " + str(monomer_mean_len))
							print("number of monomers = " + str(monomer_count))
					else:
						out_no_tandcirc_list.append(line)
						out_no_tandcirc_list.append(seq_line)
						out_no_tandcirc_list.append(end_repeat_kmer)
						if verb_debug:
							print("-----=====-----")
							print("non-tandem circle found")
							print("seqID = " + str(line))
							print("input seq        = " + str(seq_line))
							print("end repeat kmer    = " + end_repeat_kmer)
				else:
					out_no_tandcirc_list.append(line)
					out_no_tandcirc_list.append(seq_line)
					out_no_tandcirc_list.append(end_repeat_kmer)
					if verb_debug:
						print("-----=====-----")
						print("non-tandem circle found")
						print("seqID = " + str(line))
						print("input seq        = " + str(seq_line))
						print("end repeat kmer    = " + end_repeat_kmer)
		##
		return out_no_tandcirc_list, out_yes_tandcirc_list
	##
	
	if tandem_circ:
		no_tandem_circ_list, yes_tandem_circ_list = tandemcircle(circ_list)
		circ_list = no_tandem_circ_list + yes_tandem_circ_list
		if debug:
			numb_no_tandem_seq = str(int(len(no_tandem_circ_list)/3))
			numb_yes_tandem_seq = str(int(len(yes_tandem_circ_list)/3))
			print("-----=====-----")
			print("running tandem 'circle' discovery on " + numb_output_seq + " filtered 'circular' sequences")
			print("found " + numb_no_tandem_seq + " non-tandem simple or extended 'circles'")
			print("found " + numb_yes_tandem_seq + " tandem simple or extended 'circles'")
	else:
		if debug:
			print("-----=====-----")
			print("NOT running tandem 'circle' discovery on " + numb_output_seq + " filtered 'circular' sequences")
	##
	
	## rename the rnaSPAdes seqIDs to include the new seqID prefix, the new unit length, and the corrected coverage
		## also preserve the NODE ID
	
	def renamer(in_list, circ_bool):
		out_rename_fasta = []
		seq_ID_counter = 0
		digit_count = len(str(sum('>' in string for string in in_list)))
		for line_ind, line in enumerate(in_list):
			if line.startswith('>'):
				NODE_ID = '_'.join(line.split('_')[0:2])
				unit_len = int(len(in_list[line_ind+1]))
				if (unit_len >= min_len):
					start_len = int(line.split('_')[3])
					start_cov = float(line.split('_')[5])
					new_cov = round(start_cov * (start_len / unit_len),4)
					seq_ID_str = str(seq_ID_counter).zfill(digit_count)
					if circ_bool:
						new_seqID = ">" + new_seqID_prefix + "_" + NODE_ID[1:] + "_UL_" + str(unit_len) + "_AC_" + str(new_cov) + "_ID_" + seq_ID_str + "_C"
					else:
						new_seqID = ">" + new_seqID_prefix + "_" + NODE_ID[1:] + "_UL_" + str(unit_len) + "_AC_" + str(new_cov) + "_ID_" + seq_ID_str + "_L"
					seq_ID_counter += 1
					out_rename_fasta.append(new_seqID)
					out_rename_fasta.append(in_list[line_ind+1])
		##
		return out_rename_fasta
	##
	
	out_circ_fasta = renamer(circ_list, 1)
	
	## extract and re-name the "linear" contigs if specified
	
	def linextractor(in_fasta_list, in_circ_ind_list):
		out_lin_fasta_list = copy.deepcopy(in_fasta_list)
		in_circ_ind_list = sorted(in_circ_ind_list, reverse=True)
		for circ_ind in in_circ_ind_list:
			if (circ_ind < len(out_lin_fasta_list)):
				out_lin_fasta_list.pop(circ_ind)
		##
		return out_lin_fasta_list
	##
	
	if keep_linear:
		out_lin_fasta = linextractor(input_fasta_list, circ_ind_list)
		out_lin_fasta = renamer(out_lin_fasta, 0)
		if debug:
			print("-----=====-----")
			print("'linear' contigs will be parsed an saved to their own .fasta")
		##
	
	## save found circles (and linears) and report findings
	
	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##
	
	print("-----=====-----")
	
	if simple_circ:
		print("by looking for a simple k-mer repeat of length " + str(k_len))
	else:
		print("by looking for a 5' extendable k-mer repeat of seed length " + str(k_len))
	##
	
	if tandem_circ:
		print("and by further looking for tandemly repeated monomeric lengths allowing +/-" + str(round(tandem_thr,2)*100) + "% length variation")
	##
	
	print("-----=====-----")
	print("found " + str(int(len(out_circ_fasta)/2)) + " 'circles' with unit length >=" + str(min_len) + "nt")
	print("saving 'circles' to " + output_circ_fasta_name)
	saver(output_circ_fasta_name, out_circ_fasta)
	
	if keep_linear:
		print("-----=====-----")
		print(str(int(len(out_lin_fasta)/2)) + " 'linear' contigs were also kept")
		print("saving 'linears' to " + output_lin_fasta_name)
		saver(output_lin_fasta_name, out_lin_fasta)
	##
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-i', type=str, help='input rnaSPAdes derived multi .fasta e.g.: "SRR11060619_sub_rnaSPAdes.fasta"')
	parser.add_argument('-k_len', type=int, default = 73, help='k-mer length to use for detection, this could be the last k-len listed in spades.log - default = 73')
	parser.add_argument('-simple', type=int, default = 1, help='indentify "circles" with a perfect k-mer repeat between the 5` and 3` ends ("simple" = 1) or allow the 3` derived k-mer to be 5` extended when scanning from the 5` end, this is what "cirit" does ("extended" = 0) - default = 1 = "simple"')
	parser.add_argument('-tandem', type=int, default = 0, help='attempt to identify tandemly duplicated "cicles" by extracting a consistant monomeric sequence with the k-mer repeat - default = don`t look for tandem repeats (0)')
	parser.add_argument('-tandem_thr', type=float, default = 0.1, help='if looking for tandem repeats allow this much fractional length mismatch between monomers - default = 0.1x = lengths can be within 10percent of eachother')
	parser.add_argument('-min_len', type=int, default = 10, help='minimum length of identified "circles" to report - default = 10nt')
	parser.add_argument('-keep_lin', type=int, default = 0, help='write the "linear" contigs to their own output file - default = off (0)')
	parser.add_argument('-prefix', type=str, default = 'circle', help='fixed seqID prefix string to add to identified "circles" e.g. the SRR number - default = "circle"')
	parser.add_argument('-o', type=str, default = 'found_circles.fasta', help='output .fasta filename to write the identified "circles" to - default = "found_circles.fasta"')
	parser.add_argument('-oL', type=str, default = 'found_linears.fasta', help='output .fasta filename to write the identified "linears" to - default = "found_linears.fasta"')
	parser.add_argument('-debug', type=int, default = 1, help='print more detailed results - default = 1')
	parser.add_argument('-verbose', type=int, default = 0, help='print even more detailed results - default = 0')
	
	args = parser.parse_args()
	main(args.i, args.k_len, args.simple, args.tandem, args.tandem_thr, args.min_len, args.keep_lin, args.prefix, args.o, args.oL, args.debug, args.verbose)
##







































