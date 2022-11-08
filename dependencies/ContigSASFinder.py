## the purpose of this script is to take in the the parsed output from Edgar's circUCLUST (a single .fasta file containing co-clustering sequences)
## and to try and sub-divide the constituent sequences into (arbitraily ascribed) + and - sense sequences (if they are present)
## this will be achieved as follows:
	## 0) load the cluster
	## 1) for the chosen k-mer size (default = 32) extend each sequence by k-1 nucleotides on the 3' end to allow for their circularity
	## 2) for the first sequence (the centroid), count the k-mers it shares with each other sequence, also count the RC k-mers it shares
		## counting is achieved by constructing a k-mer dictionary (key = k-mer, value = count)
			## and then iterating over it for each sequence and counting the number of k-mer matches
			## IN THIS CASE, A SET WOULD HAVE PROBABLY SUFFICED BECAUSE THE VALUES AREN'T USED
		## for each other sequence the polarity grouping to the centroid is determined by whichever k-mer count is larger
		## so the centroid sequence is arbitrarily designated as "sense"
	## 3) if all the sequences turn out to be "sense", then they are all the same polarity - write to stdout, no output is saved
	## 4) if at least one "antisense" sequence is found, write two new (multi)fastas: the "sense" and the "antisense" files
		## at this point, append to each seqID a "_S" or "_A" if a sequence was determined to be "sense" or "antisense"
		## report some counting stats on outputs

import argparse
	
def main(cluster_in_name, klen, cluster_out_name_part, rename, centroid, debug):

	## declare variables

	##cluster_in_name = 'SRR11060619_sub_cluster_010.fasta'
	##klen = 32
	##cluster_out_name_part = 'SRR11060619_sub_cluster_010'
	##rename = 1
	##centroid = 0
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

	cluster_in_list = opener(cluster_in_name)

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

	cluster_in_list = singleline(cluster_in_list)

	## add 5' k-mer to 3' end

	def kmerperm(in_list):
		list_out = []
		for line_ind, line in enumerate(in_list):
			list_out.append(line)
			if (line_ind%2):
				start_kmer = line[:(klen-1)]
				new_line = line + start_kmer
				list_out[line_ind] = new_line
		return list_out
	##

	cluster_in_list_kperm = kmerperm(cluster_in_list)

	## count k-mers (sense and antisense)

	def antisense(sequence): ## define a function "antisense" that computes the reverse complement of a string
		return sequence.replace('A','t').replace('T','a').replace('G','c').replace('C','g')[::-1].upper()
	##

	def makekmerdict(in_seq):
		kmer_dict = {}
		for kmer_ind in range(len(in_seq)-klen+1):
			kmer = in_seq[kmer_ind:kmer_ind+klen]
			if (len(kmer) == klen):
				if kmer not in kmer_dict:
					kmer_dict[kmer] = 0
				kmer_dict[kmer] += 1
		return kmer_dict
	##

	centroid_seq = cluster_in_list_kperm[1]

	centroid_sense_kmer_dict = makekmerdict(centroid_seq)

	## add the centroid seqID/seq to the "sense" output list

	sense_list_out_kperm = []

	if rename:
		new_seqID = cluster_in_list_kperm[0] + "_S"
	else:
		new_seqID = cluster_in_list_kperm[0]
	##

	sense_list_out_kperm.append(new_seqID)
	sense_list_out_kperm.append(cluster_in_list_kperm[1])

	antisense_list_out_kperm = []

	## for each sequence that is not the centroid, count the sense and antisense k-mers
	## then partition the sequence to the appropriate list based on whichever k-mer dict matched the centroid more

	def matchkmercount(refernce_kmer_dict, query_kmer_dict):
		match_count = 0
		for ref_kmer in refernce_kmer_dict:
			if ref_kmer in query_kmer_dict:
				match_count += 1
		return match_count
	##

	for line_ind, line in enumerate(cluster_in_list_kperm[2:]):
		if (line_ind%2):
			member_sense_kmer_dict = makekmerdict(line)
			sense_kmer_count = matchkmercount(centroid_sense_kmer_dict, member_sense_kmer_dict)
			member_antisense_kmer_dict = makekmerdict(antisense(line))
			antisense_kmer_count = matchkmercount(centroid_sense_kmer_dict, member_antisense_kmer_dict)
			if (sense_kmer_count >= antisense_kmer_count):
				if rename:
					new_seqID = cluster_in_list_kperm[line_ind+2-1] + "_S"
				else:
					new_seqID = cluster_in_list_kperm[line_ind+2-1]
				##
				sense_list_out_kperm.append(new_seqID)
				sense_list_out_kperm.append(cluster_in_list_kperm[line_ind+2])
			else:
				if rename:
					new_seqID = cluster_in_list_kperm[line_ind+2-1] + "_A"
				else:
					new_seqID = cluster_in_list_kperm[line_ind+2-1]
				##
				antisense_list_out_kperm.append(new_seqID)
				antisense_list_out_kperm.append(cluster_in_list_kperm[line_ind+2])
			##
	##

	## report the findings (if any) and trim off the previously appended permuted k-mer

	def kmerdeperm(in_list):
		list_out = []
		for line_ind, line in enumerate(in_list):
			list_out.append(line)
			if (line_ind%2):
				new_line = line[:-(klen-1)]
				list_out[line_ind] = new_line
		return list_out
	##

	antisense_list_out = []

	if antisense_list_out_kperm:
		sense_list_out = kmerdeperm(sense_list_out_kperm)
		antisense_list_out = kmerdeperm(antisense_list_out_kperm)
		if debug:
			print(str(int(len(cluster_in_list)/2)) + " input contigs")
			print(str(int(len(sense_list_out)/2)) + " were identifed as 'sense'")
			print(str(int(len(antisense_list_out)/2)) + " were identifed as 'antisense'")
			print("saving these contigs to new sub-clusters")
	else:
		if debug:
			print("no antisense contigs were identified")
			print("no new files will be saved")
	##

	## if a sense/antisense destinction was found, save the resulting sub-clusters

	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##

	if antisense_list_out:
		if centroid:
			if debug:
				print("saving the centroid of this S-AS cluster")
			##
			centroid_out_name = cluster_out_name_part + '_centroid.fasta'
			saver(centroid_out_name, sense_list_out[:2])
		##
		sense_out_name = cluster_out_name_part + '_sense.fasta'
		saver(sense_out_name, sense_list_out)
		antisense_out_name = cluster_out_name_part + '_antisense.fasta'
		saver(antisense_out_name, antisense_list_out)
	##
	
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-inclust', type=str, help='input circUCLUST-derived (multi).fasta e.g.: "SRR11060619_sub_rnaSPAdes_cirit_cUCLUSTer_0_SL.fasta"')
	parser.add_argument('-klen', type=int, default = 32, help='k-mer length to use for counting - default = 32')
	parser.add_argument('-outstr', type=str, default = 'sub_clustered', help='fixed filename prefix string to add if sense and antisense sub-clusters are found - default = "sub_clustered"')
	parser.add_argument('-rename', type=int, default = 1, help='append an "S" or an "A" to sense and antisense seqIDs - default = 1')
	parser.add_argument('-centroid', type=int, default = 0, help='save the centroid of the S-AS cluster to its own .fasta - default = 0')
	parser.add_argument('-debug', type=int, default = 1, help='print results - default = 1')
	
	args = parser.parse_args()
	main(args.inclust, args.klen, args.outstr, args.rename, args.centroid, args.debug)
##













