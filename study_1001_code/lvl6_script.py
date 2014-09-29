# Jim Griffin
# written 9/29/14
# Example 16s network analysis script starting with qiime/ccrepe results


# imports
import pandas as pd
import numpy as np
import shared_code.shared_functions as sf


# paths to files
data_path = '../study_1001_data/s1001_lvl6/'
results_path = '../study_1001_results/s1001_lvl6/'
# QIIME data
rep_set = '/Users/jimbo/qiime/study_1001/output/dn_otu/rep_set/study_1001_seqs_rep_set.fasta'
tax_assignments = '/Users/jimbo/qiime/study_1001/output/dn_otu/uclust_assigned_taxonomy/study_1001_seqs_rep_set_tax_assignments.txt'


# Data imports (just to shorten run time while doing this iteratively - delete later
rep_seqs = pd.read_pickle(results_path+"node_dataframe")
data = pd.io.parsers.read_table(data_path+'ccrepe_results_lvl6', sep = ' ')
edges = pd.read_pickle(data_path+"all_edges")
nodes = pd.read_pickle(data_path+"all_nodes")
filtered_edges = pd.read_pickle(data_path+"filtered_edges_q05_r80")


# import QIIME data and store in a dataframe organized by:
###
# QIIME denovoID
# full taxa string
# sequence
# each inidividual subclass
###
#rep_seqs = sf.load_otu_sequences(rep_set, tax_assignments)
#rep_seqs = sf.taxa_string_split(rep_seqs, 'taxa', depth=6)
#rep_seqs.to_pickle(results_path+"node_dataframe")


# import ccrepe data and store as edge and node table. Filter edge table based on q and r values
#edges, node_table = sf.make_node_edge_tables(data)
#edges2 = sf.edge_filter(edges, .05, .8)
#edges.to_pickle(data_path+"all_edges")
#edges2.to_pickle(data_path+"filtered_edges_q05_r80")
#node_table.to_pickle(data_path+"all_nodes")



# to do:

# Next step after getting node table and edge table. Import lvl6 abundance data from QIIME directory
# check edges, frequency, abundance, abundance by sample type, interaction density, 
# inter vs intra group interaction density and interaction density vs sequence similarity

# pseudocode
#for each pair in filtered edge list:
#    edge_list[each pair][hamming distance] = sf.hamming_distance(sequence from dataframe row name matching feature1, ...)


#print len(nodes)
#print len(rep_seqs)
# How to reconcile that I have two attempts to get all of the nodes in the dataset that return very different answers? 
