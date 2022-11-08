## this script takes in as an input a list of lists (tab delim file where each sublist is a new line)
## it then uses the "connected component" script from https://www.geeksforgeeks.org/python-merge-list-with-common-elements-in-a-list-of-lists/
## to merge any sublists that have a shared element
## after merging each remaining sublist is then saved to an individual file

## the initial usecase of this script is to merge clusters that have shared members

from collections import defaultdict
import argparse

def main(input_clusters_name, output_cluster_name):
	
	##input_clusters_name = 'SRR11060619_sub_circles_extended_clusters.seqIDs'
	##output_cluster_name = 'SRR11060619_sub_cluster'
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry.split('\t')
		##
		return open_list
	##
	
	input_clusters_lists = opener(input_clusters_name)
	
	## merge function to  merge all sublist having common elements.
	
	def merge_common(lists):
		neigh = defaultdict(set)
		visited = set()
		for each in lists:
			for item in each:
				neigh[item].update(each)
		def comp(node, neigh = neigh, visited = visited, vis = visited.add):
			nodes = set([node])
			next_node = nodes.pop
			while nodes:
				node = next_node()
				vis(node)
				nodes |= neigh[node] - visited
				yield node
		for node in neigh:
			if node not in visited:
				yield sorted(comp(node))
	##
	
	merged_clusters_lists = list(merge_common(input_clusters_lists))
	merged_clusters_lists = sorted(merged_clusters_lists, key = len, reverse = True)
	
	print("-----=====-----")
	print("of " + str(len(input_clusters_lists)) + " input clusters")
	print("produced " + str(len(merged_clusters_lists)) + " merged clusters")
	print("saving")
	
	## save merged clusters
	
	def subsaver(input_name, input_list):
		counter = 0
		cluster_count = len(str(len(input_list)))
		for element in input_list:
			cluster_counter = str(counter).zfill(cluster_count)
			full_name = output_cluster_name + "_" + cluster_counter + ".seqIDs"
			name_obj = open(full_name, "w")
			if isinstance(element, list):
				for subelement in element:
					name_obj.write(subelement + "\n")
			else:
				name_obj.write(element + "\n")
			counter += 1
		##
		name_obj.close()
	##
	
	subsaver(output_cluster_name, merged_clusters_lists)
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-i', type=str, help='input multi-cluster seqID list (tab-delim members, newline clusters) e.g.: "SRR11060619_sub_circles_extended_clusters.seqIDs"')
	parser.add_argument('-o', type=str, default = 'merged_cluster', help='output prefix string for each cluster after merging - default = "merged_cluster"')
	
	args = parser.parse_args()
	main(args.i, args.o)
##































