#!/usr/bin/python

from __future__ import division
import re
import sys
import os
from ete3 import Tree, PhyloTree

#function to skip if trees do not have trees
def check_if_tree_contains_outgroups(fam_tree, outgrp_re):
	leaf_arr = get_node_leaves(fam_tree)
	flag = check_for_outgroups(leaf_arr, outgrp_re)
	return(flag)

#get leaves under any given node
def get_node_leaves(node):
	leaf_node_arr = node.get_leaves()
	leaf_arr=list()
	for leaf in leaf_node_arr:
		leaf_arr.append(leaf.name)
	
	return(leaf_arr)

#checks for one or more outgroup sequences in the list of nodes
def check_for_outgroups(leaf_arr, outgrp_re):
	outgrp_matches=list(filter(outgrp_re.search, leaf_arr))
	if(len(outgrp_matches)==0):
		return 1
	else:
		return 0

#function for reading profile file
def read_profile_file(profile_fileName):
	profile_file=open(profile_fileName, "r")
	outgrp_regex_str = ''
	ingroup_regex_str = ''
	species_dict = {}
	for line in profile_file:
		line=line.rstrip()
		linearr=re.split(r'\s+',line)
		species_id=linearr[0]
		seq_count=int(linearr[1])
		if(seq_count==0):
			outgrp_regex_str = outgrp_regex_str + "|" + species_id
		else:
			ingroup_regex_str = ingroup_regex_str + "|" + species_id
			species_dict[species_id] = seq_count

	
	outgrp_regex_str = outgrp_regex_str[1:]
	ingroup_regex_str = ingroup_regex_str[1:]
	return([outgrp_regex_str, ingroup_regex_str, species_dict])

#function to get ingroup sequence list from the tree
def get_ingroup_sequence_list(fam_tree, ingrp_re):
	leaf_arr = get_node_leaves(fam_tree)
	ingroup_matches_arr = get_ingroup_sequences(leaf_arr, ingrp_re)
	return(ingroup_matches_arr)
	
	
#function to get filter out ingroup sequences from all the leaves of the tree
def get_ingroup_sequences(leaf_arr, ingrp_re):
	ingroup_matches_arr=list(filter(ingrp_re.search, leaf_arr))
	return(ingroup_matches_arr)

# function to get ingroup sequence pairs
def get_ingroup_sequence_pairs(ingroup_matches_arr):
	ingroup_matches_arr_len = len(ingroup_matches_arr)
	ingroup_pair_arr=list()
	for ind1 in range(0,ingroup_matches_arr_len):
		for ind2 in range(ind1+1,ingroup_matches_arr_len):
			ingroup_pair_arr.append([ingroup_matches_arr[ind1],ingroup_matches_arr[ind2]])
	return(ingroup_pair_arr)

#function for checking if given ingroup pair is true/false positive
def inspect_ingroup_pairs (fam_tree, ingroup_pair_arr, outgrp_re):
	fp_pair_count=0
	tp_pair_count=0
	for ingroup_pair in ingroup_pair_arr:
		pair_mrca = fam_tree.get_common_ancestor(ingroup_pair[0], ingroup_pair[1])
		pair_mrca_leaf_arr = get_node_leaves(pair_mrca)
		outgrp_flag = check_for_outgroups(pair_mrca_leaf_arr, outgrp_re)
		if(outgrp_flag==1):
			tp_pair_count+=1
		else:
			fp_pair_count+=1
	
	precision_val = calculate_precison_value(tp_pair_count, fp_pair_count)
	return (precision_val)

#function for calculating precision value	
def calculate_precison_value(tp_pair_count, fp_pair_count):
	precision_val = tp_pair_count/(tp_pair_count+fp_pair_count)
	return(precision_val)


#processes given family file
def process_family_tree(fam_tree_fileName, outgrp_regex_str, ingroup_regex_str, species_dict):
	
	fam_tree = PhyloTree(fam_tree_fileName, format=1)

	outgrp_re = re.compile(outgrp_regex_str)
	ingrp_re = re.compile(ingroup_regex_str)


	flag = check_if_tree_contains_outgroups(fam_tree, outgrp_re)
	if(flag==1):
		return(0)

	ingroup_matches_arr = get_ingroup_sequence_list(fam_tree, ingrp_re)
	ingroup_pair_arr = get_ingroup_sequence_pairs(ingroup_matches_arr)
	precision_val = inspect_ingroup_pairs (fam_tree, ingroup_pair_arr, outgrp_re)

	print '{0} {1}'.format(os.path.splitext(os.path.basename(fam_tree_fileName))[0], precision_val)


def read_fam_tree_files_from_directory(fam_tree_dirName, profile_fileName):
	outgrp_regex_str, ingroup_regex_str, species_dict = read_profile_file(profile_fileName)
	for fam_tree_fileName in os.listdir(fam_tree_dirName):
		process_family_tree(fam_tree_dirName+"/"+fam_tree_fileName, outgrp_regex_str, ingroup_regex_str, species_dict)

fam_tree_dirName = sys.argv[1]
profile_fileName = sys.argv[2]

read_fam_tree_files_from_directory(fam_tree_dirName, profile_fileName)
