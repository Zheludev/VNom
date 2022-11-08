## written by INZ - 11/08/22
## Stanford Unversity
## provided with no acceptance of liability or promise of functionality
## version 0.1.0

import argparse
import sys
import os
global_dir = sys.path[0]
from collections import Counter
import shutil
import glob
from os.path import exists
sys.path.insert(0, global_dir+'/dependencies')
from FASTACleanUp import main as CleanUp
from CircleFinder import main as CircleFinder
from ContigSASFinder import main as SASFinder
from CommonElementsMerger import main as Merger
from RepeatResolver import main as RepeatResolver

def main(input_prefix_string,CU_max_len,CU_min_len,CF_k_len,CF_simple,CF_tandem,CF_tandem_thr,CF_min_len,circU_ID,SAS_k_len,US_vs_all,US_id,US_cov,CC_len_thr,USL_E_val,USL_id_thr,RR_overlap_fxn,RR_unit_length_fxn):
	
	## input string
	##input_prefix_string = 'SRR11060619_sub'
	
	## starting directory
	starting_dir = os.getcwd()
	
	##================================================================================================================================================
	## check that the input_prefix_string contains exactly one '_'

	char_count = Counter(input_prefix_string)
	underscore_count = char_count['_']
	if (underscore_count != 1):
		print("-----=====-----")
		print("the input prefix string needs to have exactly one underscore within it (not on the edges)")
		print("quitting")
		quit()
	##

	##================================================================================================================================================
	## length filter the input contigs and sort in decending order
	
	## CleanUp variables
	##CU_max_len = 1000
	##CU_min_len = 10
	CU_direction = 'd'
	CU_debug = 1
	
	print("-----=====-----")
	print("filtering input contigs")
	
	## CleanUp input variables
	## input fasta, max length, min length, sort length direction (d,n,a), output fasta, debug
	
	CleanUp(input_prefix_string+'.fasta', CU_max_len, CU_min_len, CU_direction, input_prefix_string+'_filt.fasta', CU_debug)
	
	##================================================================================================================================================
	## extract 'circular' contigs and store both the 'circular' and 'linear' contigs
	
	## CircleFinder variables
	##CF_k_len = 10
	##CF_simple = 0
	##CF_tandem = 1
	##CF_tandem_thr = 0.1
	##CF_min_len = 10
	CF_keep_linears = 1
	CF_debug = 1
	CF_verb_debug = 0
	
	print("-----=====-----")
	print("running CircleFinder")
	
	## CircleFinder input variables
	## input fasta, k-mer length, simple (1) or extended (0), tandem, keep linears, prefix, circle output fasta, linear output fasta, debug, verb debug
	
	CircleFinder(input_prefix_string+'_filt.fasta', CF_k_len, CF_simple, CF_tandem, CF_tandem_thr, CF_min_len, CF_keep_linears, input_prefix_string, 	input_prefix_string+'_cir.fasta', input_prefix_string+'_lin.fasta', CF_debug, CF_verb_debug)
	
	circles_found = exists(starting_dir + '/' + input_prefix_string+'_cir.fasta')
	
	if not circles_found:
		print("-----=====-----")
		print("CircleFinder did not find any 'circles' in " + input_prefix_string+'.fasta')
		print("quitting")
		quit()
	##
	
	##================================================================================================================================================
	## sort the resulting circles (if present) and running circUCLUST as a system command
	
	## circUCLUST variables
	circU_path = global_dir+'/dependencies/circuclust'
	##circU_ID = 0.7
	circU_command = circU_path + ' -cluster ' + input_prefix_string+'_cir.fasta -id ' + str(circU_ID) + ' -fastaout ' + input_prefix_string+'_cir_UC.fa 	-tsvout ' + input_prefix_string+'_cir_UC.tsv'
	
	print("-----=====-----")
	print("running circUCLUST")
	
	## input fasta, max length, min length, sort length direction (d,n,a), output fasta, debug
	
	CleanUp(input_prefix_string+'_cir.fasta', 0, 0, CU_direction, input_prefix_string+'_cir.fasta', CU_debug)
	
	print("-----=====-----")
	
	os.system(circU_command)
	
	##================================================================================================================================================
	## parse out the identified clusters with >1 members (and quit if none present)
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry
		##
		return open_list
	##
	
	circU_tsv_out_list = opener(input_prefix_string+'_cir_UC.tsv')
	
	## find 'member' contigs (demarked with a 'H\t' prefix)
		## so by definition, if a 'member' is found then at least one cluster exists with >1 contigs
	
	circU_tsv_members_list = []
	for circU_out_entry in circU_tsv_out_list:
		if circU_out_entry.startswith('H\t'):
			circU_tsv_members_list.append(circU_out_entry)
	##
	
	## quit if no >1 clusters found
	
	print("-----=====-----")
	
	if not circU_tsv_members_list:
		print("circUCLUST did not find clusters with >1 membership in " + input_prefix_string+"_cir.fasta at " + str(circU_ID) + "% ID clustering")
		print("quitting")
		quit()
	##
	
	## identify the centroid sequences of these members (the 2nd seqID in each element)
	
	circU_seqID_centroid_list = []
	for circU_member_entry in circU_tsv_members_list:
		centroid_seqID = circU_member_entry.split('\t')[2]
		if centroid_seqID not in circU_seqID_centroid_list:
			circU_seqID_centroid_list.append(centroid_seqID)
	##
	
	print("found " + str(len(circU_seqID_centroid_list)) + " >1 occupancy clusters with circUCLUST at " + str(circU_ID) + "% ID clustering")
	
	## function to make sure a loaded fasta (as a list) is one line (element) per seqID followed by one line (element) per sequence
	
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

	## function that takes in a multi.fasta as a list and a list of seqIDs (no '>'s) and creates an output multi.fasta list of those with matching seqIDs
	
	def SequenceExtractor(in_fasta_list, in_seqID_list):
		in_fasta_list = singleline(in_fasta_list)
		out_fasta_list = []
		for seqID in in_seqID_list:
			for fasta_entry_ind, fasta_entry in enumerate(in_fasta_list):
				if fasta_entry.startswith('>'):
					if (seqID == fasta_entry[1:]):
						out_fasta_list.append(fasta_entry)
						out_fasta_list.append(in_fasta_list[fasta_entry_ind+1])
		##
		return out_fasta_list
	##
	
	## generic function for saving a list as a textfile with each element on its own line and each sub-element (if present) as tab-delimitted

	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##
	
	## load 'circles'.fasta
	
	CF_circles_fasta_list = opener(input_prefix_string+'_cir.fasta')
	
	os.mkdir('0_non_singleton_clusters')

	## for each centroid, write a new directory under a new "0_non_singleton_clusters" directory in which the seqIDs will moved and extact the sequences
	
	counter = 0
	digit_count = len(str(len(circU_seqID_centroid_list)))
	members_seqID_list = []
	all_init_clusts_names_list = []
	centroid_seqID_added = False
	for centroid_seqID in circU_seqID_centroid_list:

		## pull out centroid + member seqIDs

		for circU_member_entry in circU_tsv_members_list:
			if not centroid_seqID_added:
				members_seqID_list.append(centroid_seqID)
				centroid_seqID_added = True
			##
			if (centroid_seqID == circU_member_entry.split('\t')[2]):
				member_seqID = circU_member_entry.split('\t')[1]
				members_seqID_list.append(member_seqID)
		##

		## create a directory for each cluster in '0_non_singleton_clusters'

		cluster_count_str = str(counter).zfill(digit_count)
		cluster_ID = input_prefix_string + '_CU_clust_' + cluster_count_str
		os.chdir(starting_dir+'/0_non_singleton_clusters/')
		os.mkdir(cluster_ID)

		## write the seqID list to all_init_clusts_names_list

		all_init_clusts_names_list.append(cluster_ID)

		## write the cluster's seqIDs and sequences (as .fasta) to this directory

		os.chdir(starting_dir+'/0_non_singleton_clusters/'+cluster_ID+'/')
		saver(cluster_ID+'.seqID',members_seqID_list)
		cluster_fasta_list = SequenceExtractor(CF_circles_fasta_list, members_seqID_list)
		saver(cluster_ID+'.fasta',cluster_fasta_list)
		counter += 1
		members_seqID_list = []
		centroid_seqID_added = False
	##
	
	CF_circles_fasta_list = []
	
	os.chdir(starting_dir+'/0_non_singleton_clusters/')
	saver('non_singleton_clusters_names.list',all_init_clusts_names_list)
	
	##================================================================================================================================================
	## on these identified clusters, run Sense-Antisense Contig Finder (SASFinder)
	
	## SASFinder variables
	##SAS_k_len = 32
	SAS_rename = 0
	SAS_centroid = 1
	SAS_debug = 1
	
	print("-----=====-----")
	print("running SASFinder")
	
	## SASFinder input variables
	## input cluster multi fasta, k_len, output filename prefix, rename seqIDs, save centroids seperately, debug
	
	SAS_cluster_list = []
	for init_cluster_ID in all_init_clusts_names_list:
		SAS_input_fasta = init_cluster_ID + '.fasta'
		SAS_output_string = init_cluster_ID
		os.chdir(starting_dir+'/0_non_singleton_clusters/'+init_cluster_ID+'/')
		print("-----=====-----")
		SASFinder(SAS_input_fasta, SAS_k_len, SAS_output_string, SAS_rename, SAS_centroid, SAS_debug)

		## if SASFinder identified a sense-antisense cluster, then it wrote a 'centroid.fasta' file
		## store the list of cluster directories that had a centroid output (SAS clusters)

		if exists(SAS_output_string+'_centroid.fasta'):
			SAS_cluster_list.append(init_cluster_ID)
	##
	
	os.chdir(starting_dir+'/0_non_singleton_clusters/')
	
	## quit if no SAS clusters were found
	
	if not SAS_cluster_list:
		print("-----=====-----")
		print("SASFinder did not find any clusters containing both sense and antisense sequences")
		print("quitting")
		quit()
	else:
		print("found " + str(len(SAS_cluster_list)) + " clusters containing both sense and antisense contigs")
	##
	
	## copy these identified clusters (the multi.fasta, the centroid.fasta, and the multi.seqIDs) to a new directory
	## rename the clusters to reflect the total number of SAS clusters
	## keep a log of how clusters were renamed
	
	os.chdir(starting_dir)
	os.mkdir('1_SAS_clusters')
	
	counter = 0
	digit_count = len(str(len(SAS_cluster_list)))
	cluster_rename_log_1 = []
	renamed_SAS_cluster_list = []
	for SAS_cluster in SAS_cluster_list:
		os.chdir(starting_dir+'/1_SAS_clusters/')
		SAS_cluster_count_str = str(counter).zfill(digit_count)
		counter += 1
		new_SAS_cluster_ID = input_prefix_string + '_SAS_clust_' + SAS_cluster_count_str
		os.mkdir(new_SAS_cluster_ID)
		source_path = starting_dir+'/0_non_singleton_clusters/'+SAS_cluster+'/'+SAS_cluster
		destin_path = starting_dir+'/1_SAS_clusters/'+new_SAS_cluster_ID+'/'+new_SAS_cluster_ID

		## copy over and remame the multi.fasta
		shutil.copyfile(source_path+'.fasta', destin_path+'.fasta')
		
		## copy over and rename the centroid.fasta
		shutil.copyfile(source_path+'_centroid.fasta', destin_path+'_centroid.fasta')
		
		## copy over and rename the .seqID
		shutil.copyfile(source_path+'.seqID', destin_path+'.seqID')
		
		## make a log of this renaming
		rename_log_entry = SAS_cluster + ' renamed: ' + new_SAS_cluster_ID
		cluster_rename_log_1.append(rename_log_entry)
		renamed_SAS_cluster_list.append(new_SAS_cluster_ID)
	##
	
	saver('0_to_1_renaming_log', cluster_rename_log_1)
	
	##===============================================================================================================================================
	## now build a USEARCH index of all the 'circular' and 'linear' contigs
	
	print("-----=====-----")
	print("building USEARCH index")
	
	os.chdir(starting_dir)
	
	cat_command = 'cat ' + input_prefix_string+'_cir.fasta ' + input_prefix_string+'_lin.fasta > ' + input_prefix_string+'_cat.fasta'
	
	os.system(cat_command)
	
	## CleanUp input variables
	## input fasta, max length, min length, sort length direction (d,n,a), output fasta, debug
	
	CleanUp(input_prefix_string+'_cat.fasta', 0, 0, CU_direction, input_prefix_string+'_cat.fasta', CU_debug)
	
	## USEARCH variables
	US_path = global_dir+'/dependencies/usearch'
	US_index_command = US_path + ' -makeudb_usearch ' + input_prefix_string+'_cat.fasta -output ' + input_prefix_string+'_cat.udb'
	
	print("-----=====-----")
	os.system(US_index_command)
	
	US_index_path = starting_dir + '/' + input_prefix_string+'_cat.udb'
	
	##===============================================================================================================================================
	## now run USEARCH for each SAS cluster
	## run this either in 'delta mode' (all cluster members vs index) or 'epsilon mode' (centroids only vs index)
	
	print("-----=====-----")
	print("running USEARCH vs each SAS cluster")
	
	## USEARCH variables
	##US_vs_all = 0
	##US_id = 0.7
	##US_cov = 0.8
	
	## function to append new, non-redundant seqIDs found by USEARCH_global

	def US_tsv_parse(in_tsv_output, in_seqID_list):
		for US_tsv_entry in in_tsv_output:
			seqID = US_tsv_entry.split('\t')[1]
			if seqID not in in_seqID_list:
				in_seqID_list.append(seqID)
		##
		return in_seqID_list
	##
	
	os.chdir(starting_dir+'/1_SAS_clusters/')
	
	extended_seqIDs_list = []
	extended_seqIDs_temp_list = []
	SAS_centroid_seqIDs_list = []
	for SAS_cluster in renamed_SAS_cluster_list:
		os.chdir(starting_dir+'/1_SAS_clusters/'+SAS_cluster+'/')
		US_large_tsv_out = SAS_cluster+'_USL.tsv'
		US_small_tsv_out = SAS_cluster+'_USS.tsv'

		## choose if the centroids or all members are to be used as queries
		if US_vs_all:
			print("will query all cluster members against the index")
			input_fasta = SAS_cluster+'.fasta'
		else:
			print("will query only centroids against the index")
			input_fasta = SAS_cluster+'_centroid.fasta'
		##

		print("-----=====-----")
		print("starting enrichment on " + SAS_cluster)
		US_search_command_large = US_path + ' -usearch_global ' + input_fasta + ' -db ' + US_index_path + ' -id ' + str(US_id) + ' -query_cov ' + str(US_cov) + ' -strand both -userout ' + US_large_tsv_out + ' -userfields query+target+id'
		US_search_command_small = US_path + ' -usearch_global ' + input_fasta + ' -db ' + US_index_path + ' -id ' + str(US_id) + ' -target_cov ' + str(US_cov) + ' -strand both -userout ' + US_small_tsv_out + ' -userfields query+target+id'
		
		## run USEARCH_global
		print("-----=====-----")
		os.system(US_search_command_large) ## finds contigs that are UL or >UL in length
		print("-----=====-----")
		os.system(US_search_command_small) ## finds contigs that are UL or <UL in length

		## record the seqIDs (old and those found by USEARCH)
		## want all the seqIDs per cluster to be tab-delimited element of a list

		init_seqIDs = opener(SAS_cluster+'.seqID')
		init_seqID_count = len(init_seqIDs)
		extended_seqIDs_temp_list.extend(init_seqIDs)

		large_US_tsv_out = opener(US_large_tsv_out)
		extended_seqIDs_temp_list = US_tsv_parse(large_US_tsv_out, extended_seqIDs_temp_list)

		small_US_tsv_out = opener(US_small_tsv_out)
		extended_seqIDs_temp_list = US_tsv_parse(small_US_tsv_out, extended_seqIDs_temp_list)

		ext_seqID_count = len(extended_seqIDs_temp_list)

		added_seqID_count = int(ext_seqID_count - init_seqID_count)

		extended_seqIDs_temp_list = '\t'.join(extended_seqIDs_temp_list) ## for each member, record on one line (element) with tab-delimitation
		extended_seqIDs_list.append(extended_seqIDs_temp_list) ## so each line is therefore a different cluster
		extended_seqIDs_temp_list = []

		print("found " + str(added_seqID_count) + " new sequences for " + SAS_cluster)
		
		## produce a list of centroid seqIDs from these SAS clusters

		centroid_seqID = init_seqIDs[0]
		SAS_centroid_seqIDs_list.append(centroid_seqID)
	##
	
	os.chdir(starting_dir+'/1_SAS_clusters/')
	saver('US_ext_SAS_clusters.seqIDlist', extended_seqIDs_list)
	saver('SAS_clusters_centroids.seqIDs', SAS_centroid_seqIDs_list)
	
	##================================================================================================================================================
	## for the resulting 'US_ext_SAS_clusters.seqIDlist' run Merger to merge any overlapping clusters
	
	print("-----=====-----")
	print("running Merger")
	
	os.chdir(starting_dir)
	os.mkdir('2_merged_clusters')
	source_path = starting_dir+'/1_SAS_clusters/US_ext_SAS_clusters.seqIDlist'
	destin_path = starting_dir+'/2_merged_clusters/US_ext_SAS_clusters.seqIDlist'
	shutil.copyfile(source_path, destin_path)
	
	os.chdir(starting_dir+'/2_merged_clusters/')
	
	## Merger variables
	M_out_prefix = input_prefix_string + '_merged_cluster'
	
	Merger('US_ext_SAS_clusters.seqIDlist', M_out_prefix)
	
	## find names of merged clusters
	
	merged_cluster_ID_list = []
	for filename in os.listdir(starting_dir+'/2_merged_clusters/'):
		fullpath = os.path.join(starting_dir+'/2_merged_clusters/', filename)
		if os.path.isfile(fullpath):
			if filename.endswith('.seqIDs'):
				merged_cluster_ID_list.append(filename[:-7])
	##
	
	merged_cluster_ID_list.sort()
	
	## load the catted circular and linear contigs:
	
	os.chdir(starting_dir)
	cat_contig_fasta_list = opener(input_prefix_string+'_cat.fasta')
	
	## create a directory for each merged cluster, extract the consituent sequences, and run SASFinder
	## and check if each merged cluster was indeed merged iff then choose the shorter centroid as the new centroid
	
	def CentroidCounter(centroid_seqID_list, members_seqID_list):
		centroid_count = 0
		centroid_out_list = []
		for member_seqID in members_seqID_list:
			for centroid_seqID in centroid_seqID_list:
				if (member_seqID == centroid_seqID):
					centroid_count += 1
					centroid_out_list.append(centroid_seqID)
		##
		return centroid_out_list
	##
	
	## SASFinder variables
	SAS2_k_len = SAS_k_len
	SAS2_rename = 1
	SAS2_centroid = 1
	SAS2_debug = 1
	
	for merged_cluster_ID in merged_cluster_ID_list:
		os.chdir(starting_dir+'/2_merged_clusters/')
		merged_cluster_seqID_list = opener(merged_cluster_ID+'.seqIDs')
		multi_centroid_seqIDs_list = CentroidCounter(SAS_centroid_seqIDs_list, merged_cluster_seqID_list)

		## if more than one cluster was merged:
		if (len(multi_centroid_seqIDs_list) > 1):
			print(merged_cluster_ID + " was made by merging several clusters")
			print("choosing the shortest centroid as the new centroid from the merged centroids")
		##

		## because Merger re-sorts the seqIDs, need to place the centroid seqID back on top for SASFinder to work
		## in the case where there is more than one centroid, this loop will choose the shortest to be the new centroid
		## if CU_max_len = 0, need to set shortest length really high
		if not CU_max_len:
			shortest_length = 1000000000
		else:
			shortest_length = CU_max_len
		##
		for centroid in multi_centroid_seqIDs_list:
			centroid_unit_length = int(centroid.split('_')[5])
			if (centroid_unit_length < shortest_length):
				shortest_length = centroid_unit_length
				shortest_centroid = centroid
		##
		merged_cluster_seqID_list.remove(shortest_centroid)
		merged_cluster_seqID_list.insert(0,shortest_centroid)
		##

		os.mkdir(starting_dir+'/2_merged_clusters/'+merged_cluster_ID)
		os.chdir(starting_dir+'/2_merged_clusters/'+merged_cluster_ID)
		shutil.move(starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'.seqIDs', starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'/'+	merged_cluster_ID+'.seqIDs')
		merged_fasta_list = SequenceExtractor(cat_contig_fasta_list, merged_cluster_seqID_list)
		saver(merged_cluster_ID+'.fasta',merged_fasta_list)
		
		print("-----=====-----")
		print("running SASFinder on " + merged_cluster_ID)
		SAS2_input_fasta = merged_cluster_ID+'.fasta'
		SAS2_output_string = merged_cluster_ID
		SASFinder(SAS2_input_fasta, SAS2_k_len, SAS2_output_string, SAS2_rename, SAS2_centroid, SAS2_debug)
	##
	
	cat_contig_fasta_list = []
	os.chdir(starting_dir+'/2_merged_clusters/')
	saver('merged_cluster_IDs',merged_cluster_ID_list)
	
	##================================================================================================================================================
	## lastly, the resulting clusters need to be assessed for concatemerization using MARS and USEARCH_local
	## concatemers are resolved using RepeatResolver which parses the USEARCH_local tsv output
	## this is done in a strand-aware way
	## the way SASFinder works, the centroid is defined as 'sense' so an 'antisense' centroid needs to be computed
	## first, identify any concatemeric contigs and partition them appropriately
	
	print("-----=====-----")
	print("identifying concatemeric contigs")
	
	os.chdir(starting_dir)
	os.mkdir('3_resolved_clusters')
	
	## move and rename the .fastas to '3_resolved_clusters'
	## keep a renaming record
	
	counter = 0
	digit_count = len(str(len(merged_cluster_ID_list)))
	cluster_rename_log_2 = []
	resolved_ID_cluster_list = []
	for merged_cluster_ID in merged_cluster_ID_list:
		merged_cluster_count_str = str(counter).zfill(digit_count)
		counter += 1
		new_resolved_cluster_ID = input_prefix_string+'_resolved_cluster_'+merged_cluster_count_str
		new_resolved_cluster_path = starting_dir+'/3_resolved_clusters/'+new_resolved_cluster_ID+'/'
		
		os.mkdir(new_resolved_cluster_path)

		## copy over and rename the all catted sequences
		shutil.move(starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'/'+merged_cluster_ID+'.fasta',new_resolved_cluster_path+new_resolved_cluster_ID+	'.fasta')

		## copy over and rename the antisense sequences
		shutil.move(starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'/'+merged_cluster_ID+'_antisense.fasta',new_resolved_cluster_path+	new_resolved_cluster_ID+'_antisense.fasta')

		## copy over and rename the sense sequences
		shutil.move(starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'/'+merged_cluster_ID+'_sense.fasta',new_resolved_cluster_path+	new_resolved_cluster_ID+'_sense.fasta')

		## copy over and rename the centroid sequence
		shutil.move(starting_dir+'/2_merged_clusters/'+merged_cluster_ID+'/'+merged_cluster_ID+'_centroid.fasta',new_resolved_cluster_path+	new_resolved_cluster_ID+'_centroid.fasta')

		## keep renaming log

		rename_log_entry = merged_cluster_ID + ' renamed: ' + new_resolved_cluster_ID
		cluster_rename_log_2.append(rename_log_entry)
		resolved_ID_cluster_list.append(new_resolved_cluster_ID)
	##
	
	os.chdir(starting_dir+'/3_resolved_clusters/')
	saver('2_to_3_renaming_log', cluster_rename_log_2)
	
	## go through each cluster, and by comparing to the centroid, identify any concatemeric contigs (CCs):
	## where: len(CC) > thr * len(centroid)
	## if CCs are found, write a new directory containing a multi.fasta of the CC _followed_ by the centroid
	## make sure the centroid is in the polarity of the CC (RC if needed)
	## then compute the circular permutation of the centroid vs its CC
	
	## CCFinder variables:
	##CC_len_thr = 0.1
	
	## MARS variables:
	## alphabet, input fasta, output fasta, mode, threads
	## e.g. mars -a DNA -i cat.fasta -o cp_cat.fasta -m 1 -T 24
	
	## USEARCH_local variables
	## local, query (permuted centroid), database (CC), strand polarity, e-val, id-thr, aln_out, tsv name, tsv fields
	## e.g. $IZ_py_scripts/../usearch -usearch_local pUL.fasta -db dnUL.fasta -strand plus -evalue 1e-15 -id 0.8 -alnout cp_cat.aln -userout cp_cat.tsv 	-userfields target+id+tlo+thi+trow+ql+tl
	
	USL_path = US_path + ' -usearch_local '
	##USL_E_val = '1e-15'
	##USL_id_thr = 0.8
	USL_tsv_fields = '-userfields target+id+tlo+thi+trow+ql+tl'
	
	## RepeatResolver variables
	##RR_overlap_fxn = 0.05
	##RR_unit_length_fxn = 0.1
	RR_debug = 1
	
	def ReverseComplement(sequence):
		return sequence.replace('A','t').replace('T','a').replace('G','c').replace('C','g')[::-1].upper()
	##
	
	def CCFinder(in_fasta_filename, in_centroid_filename, sense):
		## save the starting path
		initial_dir = os.getcwd()

		## excise the '.fasta' from the filename to use the remaining string for naming
		fasta_filename_start = in_fasta_filename[:-6]

		## open the .fastas needed
		cluster_fasta = opener(in_fasta_filename)
		centroid_fasta = opener(in_centroid_filename)

		## determine the upper and lower limits for considering a sequence as 'unit-length'
		centroid_upper_len_limit = (1+CC_len_thr) * int(centroid_fasta[0].split('_')[5])
		centroid_lower_len_limit = (1-CC_len_thr) * int(centroid_fasta[0].split('_')[5])

		## empty lists to load .fastas into
		CC_seqIDs = []
		UL_seqIDs = []
		sub_UL_seqIDs = []

		## determine the string and therefore polarity to use
		polarity_str = 'antisense'
		if sense:
			polarity_str = 'sense'
		##

		print("-----=====-----")

		## if working in the antisense, need to comptute the RC of the centroid for strand-aware USEARCH_local alignment
		if not sense:
			print("computing the antisense of the centroid sequence (which is sense by definition)")
			centroid_fasta[1] = ReverseComplement(centroid_fasta[1])
		##

		## scan through all cluster members and partition them based on size
		## bigger than unit-length = concatermeric contig (CC)
		## smaller than unit-length = fragment
		for entry in cluster_fasta:
			if entry.startswith('>'):
				member_len = int(entry.split('_')[5])
				if (member_len > centroid_upper_len_limit):
					CC_seqIDs.append(entry.split('>')[1]) ## this is needed bc SequenceExtractor uses seqIDs w/o '>'
				elif (member_len < centroid_lower_len_limit):
					sub_UL_seqIDs.append(entry.split('>')[1]) ## this is needed bc SequenceExtractor uses seqIDs w/o '>'
				else:
					UL_seqIDs.append(entry.split('>')[1]) ## this is needed bc SequenceExtractor uses seqIDs w/o '>'
		##

		## save sequences identified as unit-length (so omitting any >UL or <UL that were found)
		print("saving sequences identified as unit-length")
		UL_fasta = SequenceExtractor(cluster_fasta,UL_seqIDs)
		saver(fasta_filename_start + '.fasta',UL_fasta)

		## save sequences identified as sub-unit-length even prior to the RepeatResolver pipeline
		if sub_UL_seqIDs:
			print("found some initial sub unit-length sequences - saving as 'fragments'")
			sub_UL_fasta = SequenceExtractor(cluster_fasta,sub_UL_seqIDs)
			saver(fasta_filename_start + '_FR.fasta',sub_UL_fasta)
		##	

		## if CCs were found, run the MARS-USEARCH_local-RepeatResolver pipeline
		if CC_seqIDs:
			counter = 0
			digit_count = len(str(len(CC_seqIDs)))
			print("found " + str(len(CC_seqIDs)) + " concatemeric contig(s) in the " + polarity_str + " cluster of " + fasta_filename_start)
			print("creating sub-directories and populating with concatemeric contig - centroid .fasta pairs")

			## run the MARS-USEARCH_local-RepeatResolver pipeline for each pair
			for CC_contig_seqID in CC_seqIDs:
				## for each CC, create a new dir (nested by polarity) and save multi.fastas of the CC followed by the centroid
				## also save just the CC
				CC_contig_count_str = str(counter).zfill(digit_count)
				counter += 1

				## create a CC_contigs main directory to house all the sub-directories of the CC-centroid pairs
				polarity_dir_name = 'CC_contigs_' + polarity_str
				if not os.path.isdir(initial_dir+'/'+polarity_dir_name):
					os.mkdir(polarity_dir_name)
				##

				## create the CC-centroid pair directory in the main CC_contigs directory
				CC_dir_name = 'CCCpair_'+polarity_str+'_'+CC_contig_count_str
				os.mkdir(polarity_dir_name + '/' + CC_dir_name)

				## create the list of the di.fasta of the CC + the centroid contig
				CC_centroid_pair_fasta = SequenceExtractor(cluster_fasta,[CC_contig_seqID]) + centroid_fasta ## SequenceExtractor takes in lists
				
				## save this CC + centroid contig list in the new sub-directory
				current_dir = os.getcwd()
				CCCpair_path = current_dir + '/' + polarity_dir_name +'/' + CC_dir_name + '/'
				saver(CCCpair_path + CC_dir_name + '.fasta',CC_centroid_pair_fasta)

				## also save just the concatermeric contig as its own file
				saver(CCCpair_path + CC_dir_name + '_CC.fasta',CC_centroid_pair_fasta[:2])

				## now run MARS on this CC-centroid pair (thus permuting the centroid vs the CC)
				print("-----=====-----")
				print("running MARS on " + CC_dir_name)
				print("-----=====-----")
				os.chdir(CCCpair_path)
				MARS_command = global_dir+'/dependencies/MARS/mars -a DNA -i ' + CC_dir_name + '.fasta -o ' + CC_dir_name + '_CP.fasta -m 1 -T 24'
				os.system(MARS_command)

				## MARS will fail on highly degenerate circles (-P too high error)
				## in which case, skip these sequences as this is a bad sign

				MARS_output = exists(CCCpair_path + CC_dir_name + '_CP.fasta')

				if MARS_output:
					## now extract the circularly permuted centroid, rename its seqID to reflect the permutation, and save a new .fasta
					CP_centroid_fasta = opener(CC_dir_name + '_CP.fasta')[2:]
					CP_centroid_fasta[0] = CP_centroid_fasta[0] + '_CP'
					saver(CC_dir_name + '_CP_centroid.fasta',CP_centroid_fasta)
	
					## now compute the USEARCH local alignment of the permuted centroid vs the CC
					print("-----=====-----")
					print("running USEARCH_local on " + CC_dir_name)
					print("-----=====-----")
					USL_command = USL_path + CC_dir_name + '_CP_centroid.fasta -db ' + CC_dir_name + '_CC.fasta -strand plus -evalue ' + USL_E_val + ' -id ' +	 str(USL_id_thr) + ' -alnout ' + CC_dir_name + '_USL.aln -userout ' + CC_dir_name + '_USL.tsv ' + USL_tsv_fields
					os.system(USL_command)

					## USEARCH_local will output a .tsv even if it's empty, so need to check if results were found
					USL_output_size = os.stat(CCCpair_path + CC_dir_name + '_USL.tsv').st_size

					if USL_output_size:
						## now parse out the USEARCH_local result with RepeatResolver
						print("-----=====-----")
						print("running RepeatResolver on " + CC_dir_name)
						RepeatResolver(CC_dir_name + '_USL.tsv', fasta_filename_start, CC_contig_count_str, RR_overlap_fxn, RR_unit_length_fxn, RR_debug)
		
						## now if any files were written, they'd contain the 'fasta_filename_start' string in their name
						## use glob to copy the files up two directories
						for RR_out_file in glob.glob(CCCpair_path+"/"+fasta_filename_start+"*"):
							shutil.copy(RR_out_file, initial_dir)
						##
					else:
						print("-----=====-----")
						print("USEARCH_local appears to have written an empty .tsv output, indicating no hits being found")
						print("skipping these sequences as this is a bad sign")
					##
				else:
					print("-----=====-----")
					print("MARS probably experienced a 'too large -P error' - this is likely due to high degeneracy")
					print("skipping these sequences as this is a bad sign")
				##

				## return to starting directory
				os.chdir(initial_dir)
		else:
			print("did not find concatemeric contigs in the " + polarity_str + " cluster of " + fasta_filename_start)
	##
	
	for resolved_cluster_ID in resolved_ID_cluster_list:
		os.chdir(starting_dir+'/3_resolved_clusters/'+resolved_cluster_ID+'/')
		print("-----=====-----")
		print("assessing the sense hits of " + resolved_cluster_ID + " for concatemeric and sub unit-length contigs")
		CCFinder(resolved_cluster_ID+'_sense.fasta',resolved_cluster_ID+'_centroid.fasta',1)
		print("-----=====-----")
		print("assessing the antisense hits of " + resolved_cluster_ID + " for concatemeric and sub unit-length contigs")
		CCFinder(resolved_cluster_ID+'_antisense.fasta',resolved_cluster_ID+'_centroid.fasta',0)
	##
	
	##================================================================================================================================================
	##================================================================================================================================================
	## finally, cat any newly identified unit length contigs to their appropriate multi-fasta and copy these ULs + fragments to an output directory
	
	os.chdir(starting_dir)
	os.mkdir('4_final_clusters')
	
	## create the final directory and sub-directories for sense and antisense sequences
	
	os.chdir(starting_dir+'/4_final_clusters/')
	final_cluster_ID_list = []
	cluster_rename_log_3 = []
	for resolved_cluster_ID in resolved_ID_cluster_list:

		## extract the final cluster ID number from the previous directory names
		final_cluster_name = input_prefix_string + "_cluster_" + resolved_cluster_ID.split('_')[-1]

		## make nested directories
		os.mkdir(final_cluster_name)
		os.mkdir(final_cluster_name+'/sense')
		os.mkdir(final_cluster_name+'/antisense')

		## keep a list of the new, final directory names
		final_cluster_ID_list.append(final_cluster_name)

		## create a renaming log
		rename_log_entry = resolved_cluster_ID + ' renamed: ' + final_cluster_name
		cluster_rename_log_3.append(rename_log_entry)
	##
	
	os.chdir(starting_dir+'/4_final_clusters/')
	saver('3_to_4_renaming_log', cluster_rename_log_3)
	
	## ok copy over the appropriate (ULs and fragments) .fasta files to the final directory
	## funtion just to not duplicate the code for sense and antisense
	
	def tofinalcopy(polarity):
		## uses glob.glob to look for files matching a regex string
		wc_filepath = starting_dir+'/3_resolved_clusters/'+resolved_cluster_ID+'/'+resolved_cluster_ID+'_'+polarity+'*'
		for wc_fasta in glob.glob(wc_filepath):
			destin_path = starting_dir+'/4_final_clusters/'+final_cluster_ID_list[resolved_cluster_ID_ind]+'/'+polarity+'/'
			shutil.copy(wc_fasta, destin_path)
	##
	
	for resolved_cluster_ID_ind, resolved_cluster_ID in enumerate(resolved_ID_cluster_list):
		os.chdir(starting_dir+'/3_resolved_clusters/'+resolved_cluster_ID+'/')
		tofinalcopy('sense')
		tofinalcopy('antisense')
	##
	
	## cat and rename the final .fasta files
	## funtion just to not duplicate the code for sense and antisense
	
	def finalrename(polarity):
		## uses glob.glob to look for files matching a regex string

		## cat all unit-length sequences, rename to the final naming convention, and then remove the old copies
		os.chdir(starting_dir+'/4_final_clusters/'+final_cluster_ID+'/'+polarity+'/')
		if glob.glob(starting_dir+'/4_final_clusters/'+final_cluster_ID+'/'+polarity+'/*UL.fasta'):
			os.system('cat '+resolved_ID_cluster_list[final_cluster_ID_ind]+'_'+polarity+'.fasta *UL.fasta > '+final_cluster_ID+'_'+polarity+'.fasta')
			os.system('rm '+resolved_ID_cluster_list[final_cluster_ID_ind]+'_'+polarity+'.fasta')
			os.system('rm *UL.fasta')
		else:
			os.system('mv '+resolved_ID_cluster_list[final_cluster_ID_ind]+'_'+polarity+'.fasta '+final_cluster_ID+'_'+polarity+'.fasta')
		##

		## rename all fragment sequences to the final naming convention
		fragment_filepath = starting_dir+'/4_final_clusters/'+final_cluster_ID+'/'+polarity+'/*_FR.fasta'
		fragment_counter = 0
		for fragment in glob.glob(fragment_filepath):
			fragment_filename = final_cluster_ID+'_'+polarity+'_fragment_'+str(fragment_counter)+'.fasta'
			fragment_counter += 1
			os.system('mv '+fragment+' '+fragment_filename)
		##

		## give some metrics of how many sequences were saved in the end
		final_cluster_list = opener(final_cluster_ID+'_'+polarity+'.fasta')
		seq_count = int(len(final_cluster_list)/2)
		print("-----=====-----")
		print("found: " + str(seq_count) + " unit-length sequences in "+ final_cluster_ID+"_"+polarity+".fasta")
	##
	
	os.chdir(starting_dir+'/4_final_clusters/')
	for final_cluster_ID_ind, final_cluster_ID in enumerate(final_cluster_ID_list):
		finalrename('sense')
		finalrename('antisense')
	##
	
	## ===============================================================================================================================================
## 
## ok here are all the variables that make sense for the user to have control over:

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', type=str, help='input prefix string of the file to be processed - will be used in naming all output files e.g.: from the input file "SRR11060619_sub.fasta" declare "SRR11060619_sub"')
	parser.add_argument('-max', type=int, default = 1000, help='maximum (inclusive) sequence length to be considered from input - default = 1000nt')
	parser.add_argument('-min', type=int, default = 10, help='minimum (inclusive) sequence length to be considered from input - default = 10nt')
	parser.add_argument('-CF_k', type=int, default = 10, help='initial k-mer length to be used in CircleFinder - default = 10')
	parser.add_argument('-CF_simple', type=int, default = 0, help='bool on if "simple" (1) or "extended" (0) CircleFinder should be run - default = 0 ("extended"). "simple" looks for a trivial duplication of the 3prime end at the 5prime end of length -CF_k. "extended" scans from the 5prime end inwards for a match to the 3prime end of length -CF_k and then "extends" the match to the 5prime start and looks for that "extended" k-mer match at the 3prime end - in effect "extended" is the "simple" method for when the user does not know what length of -CF_k with which to start')
	parser.add_argument('-CF_tandem', type=int, default = 1, help='bool on if "tandem" (1) or not (0) CircleFinder should be run - default = 1 ("tandem"). "tandem" attempts to identify any tandem duplications within a contig that are demarked by the indentified ("simple" or "extended") 3prime k-mer, if the lengths of the duplications make sense, each monomer of the duplication is saved as its own entry in the output .fasta')
	parser.add_argument('-CF_tandem_thr', type=float, default = 0.1, help='used if "tandem" mode is selected. states the +/- length threshold permissable for a putative tandem duplication to be considered correct - default = 0.1 = tandem duplications need to be within +/- 10 percent length of eachtother')
	parser.add_argument('-CF_min', type=int, default = 10, help='minimum (inclusive) sequence length of identified "circles" to be saved to output - default = 10nt')
	parser.add_argument('-CU_ID', type=float, default = 0.7, help='fractional sequence identity threshold to be used during circUCLUST - default = 0.7 = clusters formed with a maximal pairwise identity to the centroid at 70 percent')
	parser.add_argument('-SAS_k', type=int, default = 32, help='k-mer length to be used during the Sense-AntisenseFinder cluster filtering - default = 32')
	parser.add_argument('-USG_vs_all', type=int, default = 0, help='bool on if the cluster enrichment should be based on USEARCH_global alignments to all (1) the cluster members or just the centroids (0) - default = 0. aligning to all might increase sensitivity but might also lead to erroneous merging of distantly related clusters')
	parser.add_argument('-USG_ID', type=float, default = 0.7, help='fractional sequence identity threshold to be used during USEARCH_global cluster enrichment - default = 0.7 = hits to clusters will be considered if they share at least 70 percent sequence identity with the cluster member they are aligned to')
	parser.add_argument('-USG_cov', type=float, default = 0.8, help='fractional sequence coverage threshold to be used during USEARCH_global cluster enrichment - default = 0.8 = hits to clusters will be considered if they share at least 80 percent sequence coverage with the cluster member they are aligned to')
	parser.add_argument('-CC_len_thr', type=float, default = 0.1, help='states the length threshold permissable for a concatemeric contig identification as well as sub-unit-length contig identification - default = 0.1 = unit length contigs are thus defined as within +/- 10 percent the centroid sequence length. shorter are the sub-unit-length contigs. longer are the concatemeric contigs parsed by RepeatResolver')
	parser.add_argument('-USL_Eval', type=str, default = '1e-15', help='standard form notation of the E-value threshold imposed on USEARCH_local alignment during concatemeric contig resolution vs the centroid contig - default = 1e-15 = hits to the concatemeric contig will only be considered if their alignment to the centroid scores an E-value at most 1e-15')
	parser.add_argument('-USL_ID', type=float, default = 0.8, help='fractional sequence identity threshold to be used during USEARCH_local alignment during concatemeric contig resolution vs the centroid contig - default = 0.8 = hits to the concatemeric contig will only be considered if they share at least 80 percent sequence identity with the centroid')
	parser.add_argument('-RR_overlap_thr', type=float, default = 0.05, help='fractional sequence length of overlaps or spaces permissable between hits following USEARCH_local alignment to be parsed by RepeatResolver - default = 0.05 = hits need to overlap to be spaced apart by no more than 5 percent of the centroid sequence length')
	parser.add_argument('-RR_UL_thr', type=float, default = 0.1, help='fractional sequence length threshold for putative unit-length hits to be considered following USEARCH_local alignment to be parsed by RepeatResolver - default = 0.1 = RepeatResolver will consider a hit unit-length if it is within +/- 10 percent of the length of the centroid contig')	
	args = parser.parse_args()
	main(args.i, args.max, args.min, args.CF_k, args.CF_simple, args.CF_tandem, args.CF_tandem_thr, args.CF_min, args.CU_ID, args.SAS_k, args.USG_vs_all, args.USG_ID, args.USG_cov, args.CC_len_thr, args.USL_Eval, args.USL_ID, args.RR_overlap_thr, args.RR_UL_thr)
##



























