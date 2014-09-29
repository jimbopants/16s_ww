# Data formatting functions for generating and manipulating dataframes, getting OTU names from taxa ID strings etc.


# add all the imports I need 
import pandas as pd
import numpy as np
import string



def taxa_string_split(df, name_column, sep=';', depth=5):
    """ df[name_column] -> df[kingdom, phylum, class...] to desired depth
    This assumes you have a dataframe with one (or more) name columns. 
    By default it splits a taxa string in the QIIME format:
    "k__Bacteria; p__Verrucomicrobia; c__Opitutae; o__Opitutales; f__Opitutaceae;
    g__Opitutus; s__"
    into individual columns to make comparisons at many different phylogenetic levels possible.
    """
    taxa_cat = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    for i in range(depth):
        splits = df[name_column].str.split(';').str.get(i)
        df[taxa_cat[i]] = splits

    return df



def dict_key_from_value(a_dict,a_value):
    """ returns first key that matches a dict value"""
    for key,value in a_dict.iteritems():
        if value == a_value:
            return key


def load_otu_sequences(rep_set, tax_assignments):
    """ Load files for sequences and names into a dataframe
    Returns: pd.rep_seqs with columns = ID, full taxa string, and sequence 
    """
    f = open(rep_set, 'rb')
    sequences = f.read()
    sequences = sequences.split('\n')
    f.close()
    # Names
    g = open(tax_assignments,'rb')
    seq_names = g.read()
    seq_names = seq_names.split('\n')
    g.close()

    # Extract sequence data I need
    # Takes just the "denovo#" section of the sequence names (even rows)
    seq_IDs = [ s.split(' ')[0][1:] for s in sequences[::2]]
    seq_strings = [s for s in sequences[1::2]]
    seq_IDs = seq_IDs[:-1]
    sequences = dict(zip(seq_IDs, seq_strings))  # Dict of format{ denovo#: 'AAAAT...'}
    
    # split at tabs and remove any rows that do not conform to data pattern
    seq_name_split = [s.split('\t') for s in seq_names]
    seq_name_split = [x for x in seq_name_split if len(x) >1]    
    # Get phyla strings and strip species section. 
    phyla = [ row[1] for row in seq_name_split ]
    phyla = [string.split('; s')[0] for string in phyla ]

    # Create dictionary of format {denovo# : 'k__p__c__...'} 
    names = [ row[0] for row in seq_name_split ]
    phyla_OTU_IDs = dict(zip(names, phyla))

    phyla_set = set(phyla_OTU_IDs.values())

    # Invert dictionary to only have 1 denovo assigned OTU per database aligned OTU
    shortened_dict = {}
    for taxa in phyla_set:
        if taxa not in shortened_dict:
            shortened_dict[taxa] = dict_key_from_value(phyla_OTU_IDs,taxa)
    
    rep_seqs = pd.DataFrame({'taxa' : shortened_dict.keys(), 'seq_IDs' : shortened_dict.values() })
    
    for key, val in sequences.iteritems():
        rep_seqs.loc[rep_seqs.seq_IDs==key, 'sequence'] = val

    return rep_seqs 


def make_node_edge_tables(ccrepe_results):
    # Creates and returns edge and node dataframes from ccrepe result table

    ccrepe2 = ccrepe_results.dropna() # eliminate NAs in edge table
    
    # initialize a dataframe and add rows for each unique taxon in dataset.
    node_table = pd.DataFrame( columns = ['taxa'])
    node_set = set()
    for taxa in ccrepe2.iterrows():
        for i in range(2):
            if taxa[1][i] not in node_set:
                node_set.add(taxa[1][i])
                node_series = pd.Series(taxa[1][i], index= ['taxa'])
                node_table = node_table.append(node_series, ignore_index = True)
    
    return ccrepe2, node_table


def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return (sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)))/float(len(s1))


def edge_filter(edge_table, Q_val=.05, R_val=.8):
    # Filter edges based on Q and R
    # Currently ignores negative associations which is really stupid. 
    edges = edge_table[edge_table['q.value'] <.05]
    edges = edges[edges['sim.score']>.80]
    return edges








# broken stuff below here..


# There is a function to filter dataframe below but I think it can be waay
# shorter. I don't explicitly name nodes but I can use the whole name for now I
# think. 


def retarded_version():
    # Initialize Node table:
    node_table = pd.DataFrame( columns = ['taxa', "av_abundance","edges"])
    node_list = []
    
    
    
    index_series = ['taxa', 'name', "kingdom", "phylum", "class", "order", "family", "av_abundance","edges"]
    
    # Make simplified names and taxa columns
    for match in data2.iterrows():
        row_info = []       # Intiialize rows from df entries
        row_info_2 = []
        # Get and split match info
        str_match_1, str_match_2 = str(match[1][0]), str(match[1][1]) 
        split_match_1, split_match_2 =  str_match_1.split(";"), str_match_2.split(";")
        
        # Create new row entry for dataframe filling taxa columns as well as
        # name with lowest level of classification
        for item in split_match_1:
            if len(item)>3:
                row_info.append(item[3:])
                name_1 = item
            if len(item)<4:
                row_info.append("unassigned")
        
        # Repeat for second node in assocation
        for item in split_match_2:
            if len(item)>3:
                row_info_2.append(item[3:])
                name_2 = item
            elif len(item)<4:
                row_info_2.append("unassigned")
    
        # If name is not present in name list, append it, then initialize abundance and #
        # of edges with 0's in row for now. 
        if name_1 not in node_list:
            node_list.append(name_1)
            row_info.append(0)
            row_info.append(0)
            row_info = [match[1][0]] + [name_1] + row_info
            series_row = pd.Series(row_info, index = index_series)
            node_table = node_table.append( series_row, ignore_index=True)
        if name_2 not in node_list:
            node_list.append(name_2)
            row_info_2.append(0)
            row_info_2.append(0)
            row_info_2 = [match[1][1]]+ [name_2] + row_info_2
            series_row = pd.Series(row_info_2, index = index_series)
            node_table = node_table.append( series_row, ignore_index=True)
    
    # Filter edges based on Q and R
    edge_table = data2[data2['q.value'] <.05]
    edge_table_pos = edge_table[edge_table['sim.score']>.80]
    print edge_table_pos["feature1"]
    print node_table["taxa"]
    # Convert edge_table names to match with node_table names.
    # This code is now almost identically replicated 3 times in this file, in
    # the future, don't write code like this. 
    for edge in edge_table_pos.iterrows():
        node_1,node_2 = str(edge[1][0]), str(edge[1][1])
        split_node_1, split_node_2 = node_1.split(";"), node_2.split(";")
    
        for item in reversed(split_node_1):
            if len(item)>3:
                name_1 = item
            break
        for item in reversed(split_node_2):
            if len(item)>3:
                name_2 = item
            break
        
        #edge_table_pos["feature1"][edge_table_pos["feature1"] == edge[1][0]] = name_1
        #edge_table_pos["feature2"][edge_table_pos["feature2"] == edge[1][1]] = name_2

    for nodes in edge_table_pos.iterrows():
        node_table["edges"][node_table["taxa"] == nodes[1][0]] +=1

    store_df(node_table, path, "s1001_nodes_2")
    store_df(edge_table, path, "s1001_edges_2")
    return

