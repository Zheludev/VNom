## this is a quick script that takes as input a (multi).fasta file and:
	## 1) - "single-lines" it (one line per seqID and seq)
	## 2) - length filters (upper and lower bounds allowed)
	## 3) - sorts the file based on length (both directions)

## load libraries

import argparse

def main(fasta_in_name, upper_len, lower_len, sort_dir, fasta_out_name, debug):

	## set variables

	##fasta_in_name = 'SRR11060619_sub.fasta'
	##upper_len = 1000
	##lower_len = 0
	##sort_dir = 'd'
	##fasta_out_name = 'SRR11060619_sub_clean.fasta'
	##debug = 1

	## open file

	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry
		##
		return open_list
	##

	fasta_in_list = opener(fasta_in_name)

	## make sure input sequence is singleline

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

	fasta_in_list = singleline(fasta_in_list)

	## filter sequences by length

	def lengthfilter(in_fasta_list, upper, lower):
		out_fasta_list = []
		for line_ind, line in enumerate(in_fasta_list):
			if line.startswith('>'):
				seq = in_fasta_list[line_ind+1]
				seq_len = len(seq)
				if (upper > 0) and (lower > 0):
					if (seq_len <= upper) and (seq_len >= lower):
						out_fasta_list.append(line)
						out_fasta_list.append(seq)
					##
				elif (upper > 0):
					if (seq_len <= upper):
						out_fasta_list.append(line)
						out_fasta_list.append(seq)
					##
				elif (lower > 0):
					if (seq_len >= lower):
						out_fasta_list.append(line)
						out_fasta_list.append(seq)
					##
				elif (upper == 0) and (lower == 0):
					out_fasta_list.append(line)
					out_fasta_list.append(seq)
				##
		##
		return out_fasta_list
	##

	fasta_filt_list = lengthfilter(fasta_in_list, upper_len, lower_len)

	## sort sequences by length

	def lengthsorter(in_fasta_list, direction):
		def LETdeconvolver(in_LE_list, direction):
			out_fasta_list = []
			if (direction == 'd'):
				in_LE_list = sorted(in_LE_list, key=lambda x: x[0], reverse=True)
			elif (direction == 'a'):
				in_LE_list = sorted(in_LE_list, key=lambda x: x[0], reverse=False)
			else:
				print("error, sort direction can either be 'n' (no sort - default), 'd' (descending sort), or 'a' (ascending sort)")
				print("quitting")
				quit()
			##
			for entry in in_LE_list:
				seqID = entry[1]
				seq = entry[2]
				out_fasta_list.append(seqID)
				out_fasta_list.append(seq)
			##
			return out_fasta_list
		##
		out_fasta_list = []
		temp_list = []
		length_entry_list = []
		if (direction == 'n'):
			out_fasta_list.extend(in_fasta_list)
		else:
			for line_ind, line in enumerate(in_fasta_list):
				if line.startswith('>'):
					seq = in_fasta_list[line_ind+1]
					seq_len = len(seq)
					temp_list.append(seq_len)
					temp_list.append(line)
					temp_list.append(seq)
					length_entry_list.append(temp_list)
					temp_list = []
			##
			out_fasta_list = LETdeconvolver(length_entry_list, direction)
		##
		return out_fasta_list
	##

	fasta_sort_list = lengthsorter(fasta_filt_list, sort_dir)

	## debug

	if debug:
		input_seq_count =  str(int((len(fasta_in_list)/2)))
		output_seq_count = str(int((len(fasta_sort_list)/2)))
		print("-----=====-----")
		print("input number of sequences = " + input_seq_count)

		if (upper_len > 0) and (lower_len > 0):
			print("keeping sequences above and including " + str(lower_len) + "nt")
			print("keeping sequences up to and including " + str(upper_len) + "nt")
		elif (lower_len > 0):
			print("keeping sequences above and including " + str(lower_len) + "nt")
		elif (upper_len > 0):
			print("keeping sequences up to and including " + str(upper_len) + "nt")
		elif (upper_len == 0) and (lower_len == 0):
			print("no length-based filtering performed")
		##

		if (sort_dir == 'n'):
			print("no length-based sorting permformed")
		elif (sort_dir == 'a'):
			print("sequences length-based sorted in ascending order")
		elif (sort_dir == 'd'):
			print("sequences length-based sorted in descending order")
		##

		print("output number of sequences = " + output_seq_count)
	##


	## save fasta

	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##

	saver(fasta_out_name, fasta_sort_list)
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-i', type=str, help='input (multi).fasta name e.g.: "SRR11060619_sub.fasta"')
	parser.add_argument('-upper', type=int, default = 0, help='OPTIONAL: upper bound (inclusive) sequence length limit - default = 0 (off)')
	parser.add_argument('-lower', type=int, default = 0, help='OPTIONAL: lower bound (inclusive) sequence length limit - default = 0 (off)')
	parser.add_argument('-sort', type=str, default = 'n', help='OPTIONAL: sort sequences by length - "n" no sort (default), "a" ascending, "d" descending')
	parser.add_argument('-o', type=str, default = 'cleaned.fasta', help='output (multi).fasta - default = "cleaned.fasta"')
	parser.add_argument('-debug', type=int, default = 1, help='print results - default = 1')
	
	args = parser.parse_args()
	main(args.i, args.upper, args.lower, args.sort, args.o, args.debug)
##





























