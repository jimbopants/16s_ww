"""
Written by Jim Griffin
8/14/14


Goals outlined below:

Network analysis continued:
Now that I have an association network for soil genera from the
Wilson dataset using ccrepe and have attempted to visualize this network 
with gephi, next objective is to apply basic network analysis methods to this network.

Some general questions that I hope to answer:

What is the overall network structure? Does the connectivity follow a power-law? 
Is this network a "small-world" network?
data requirements: Degree distribution, mean path length, diameter, clustering coefficient...
modularity index

Are there "hubs" within this network? Do hubs behave differently than non-hubs?
Data requirements: nodes sorted by degree, some cutoff criteria. Then
identification of within and between module hubs.  Degree vs abundance chart. 

Do phylogenetically more related groups interact differently than random groups?
Data requirements: 


"Specialists" vs generalists: are organisms that are highly abundant in one ecosystem but not present in others
different (statistically, over all node/edges) than organisms that exist in most environments?"""


##### Structure/functions
# 1. Imports
# 2. Data import
# 3. 
# 4. Generating Degree distribution plot
# 5. 


# Imports
from shared_code import association_networks
import networkx as nx
import numpy as np
import pandas as pd
from collections import Counter
from pylab import *
import math
from scipy.stats import linregress



# Useful & completed functions:
def open_saved(path="/Users/jimbo/Desktop/working/"
, results="study_1001_results/"
, data_folder="study_1001_data/"
):
    # Retrieve saved data for project. 
    
    ######     Returns:
    # path: the directory my data is stored in
    # results, raw_data, subdirectories
    # data: dataframe of correlation scores, q values, r values etc.
    # edge_list: dataframe of OTU1, OTU2, R and Q for each association
    # raw_data: contains the sample X OTU abundances for the entire dataset +
    # simplified name for OTU
    # node_list: df containing the shortened taxa name, information for that
    # taxa, number of edges, average abundance and more. 
    # meta_data: stores raw data from study 1001 mapping file.
    # connected: Nodes with >0 edges

    data = pd.io.parsers.read_csv(path+data_folder+"L5_full_table.csv", header=0, index_col=0)
    full_data = pd.io.parsers.read_table(path+data_folder+"otu_table_L5.txt",
        sep="\t")
    raw_data = pd.read_pickle(path+results+"raw_data_updated")
    edge_list = pd.read_pickle(path+results+"s1001_edges_2")
    node_list = pd.read_pickle(path+results+"s1001_nodes_2")
    connected = pd.read_pickle(path+results+"connected_nodes_df")
    meta_data = pd.io.parsers.read_table(path+data_folder+"s1001_map_file.txt", sep="\t")
    filtered_edges = pd.read_pickle(path+results+"filtered_edges")
    return (path, data, raw_data, edge_list, node_list, meta_data, connected,
full_data, filtered_edges)


def store_df(df, path, name):
    # I Would like to store dataframes as objects outside of python and
    # automatically update them. 
    # This should be deleted, it does absolutely nothing and clutters code.
    df.to_pickle(path+"study_1001_df/"+name)
    return


def get_otu_names(df, new_column):
    # Need to take a raw EMP dataframe and add simple OTU/taxa names to easily
    # match other edge/node tables.
    df[new_column] = df["Taxon"]                    # add name column
    for match in df.iterrows():                 # iterate through rows
        str_match= str(match[1][0])             # get taxa column
        split_match =  str_match.split(";")     # split into Kingdom, Phylum ...
        df.loc[match[0],new_column] = split_match[2]
        #for item in reversed(split_match):      # reverse list of taxa class
        #    if len(item)>3:                     # take first non-blank entry
        #        df.loc[match[0],"name"] = item  #Use label based indexing to set name
        #        break
    return df


# Class Interaction Patterns
def get_phyla_level(otu_string, level):
    # Splits a k;p;c;o... string and returns the desired taxa level.
    split = otu_string.split(";")
    return split[level]


def intra_class_connections(edge_list):
    # Returns a dictionary of class:# pairs that correspond to intra-class
    # connections
    phyla = {}
    for edge in edge_list.iterrows():
        if (get_phyla_level(edge[1]['feature1'],2) ==
        get_phyla_level(edge[1]['feature2'],2)):
            clasz = get_phyla_level(edge[1]['feature1'],2)
            print edge[1]['feature1'], edge[1]['feature2']
            if clasz in phyla:
                phyla[clasz]+=1
            else:
                phyla[clasz]=1
    return phyla


def calculate_class_level_connections(connected):
    # Take all of the nodes in my "connected" graph, divide into classes, calculate max
    # nodes, total nodes, and max intra-class nodes
    # Then I need to store everything as name: max, actual, max-intra, 

    edge_groups = connected.groupby("class")    # Group node list into classes
    for name, group in edge_groups:
        
        if len(group)>1:                        # Calculates actual edges
            print name#, group, len(group)
            #print group["edges"].sum()
            #print name, len(group)
            max_edges = 0
            i=1
            nodes = float(len(group)-1)
            max_intra_class = (float(len(group)) * (float(len(group))-1.0))/2.0 
            # Calculates max intra-class edges
            #print max_intra_class
            while nodes>0:
                max_edges += float(len(connected)-i)
                nodes = nodes -1
            i += 1
            print "max edges for", name, max_edges     #Calculates max number of edges
    # Need to caclulate the max intra-class edges



### Import all my data:
(path, data, raw_data, edge_list, node_list, meta_data, connected, full_data,
filtered_edges) = open_saved()



print intra_class_connections(filtered_edges)
print node_list["taxa"]

#calculate_class_level_connections(connected)

        

# Below is a column filter for pandas dataframes. I can split into
# different groups, calculate abundance in different groups etc,
# For example: Group samples by environment type
types =  meta_data["#SampleID"].groupby(meta_data["COMMON_NAME"])
for type_thing in types:
    for sample in type_thing[1]:
        for column in raw_data.columns:
            if sample == column:
                pass
                #print column, type_thing


# With column filters I can start considering generalists/specialists
def organize_raw_data(df,path):
    # This function takes a raw EMP dataframe, adds a simplified naming
    # convention, and adds metadata about each OTU. Is it present in all
    # environments or a subset? Is it highly abundant? How many total samples
    # does it appear in?
    df = get_otu_names(df) # Add names
    store_df(df,path,"raw_data_updated") # store resulting df
    return


#organize_raw_data(raw_data,path)







# These are completed but not very well written. 
def make_node_edge_tables(data, path):
    # Creating: edge/node Dataframes from otu table, goal is to get a dataframe
    # I can refer to for a while with all phylogenetic classifications and
    # network properties specified. 
    # This is the longest code to do this little ever. ROFL
    print data
    data2 = data.dropna() # eliminate NAs

    # Initialize Node table:
    node_table = pd.DataFrame( columns = ['taxa','name', "kingdom", "phylum", "class", "order", "family", "av_abundance","edges"])
    index_series = ['taxa', 'name', "kingdom", "phylum", "class", "order", "family", "av_abundance","edges"]
    node_list = []
    
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


def plot_abundance(node_table):
    # I can edit this to plot more diverse things as I need. 
    # It's good to keep this separate from what I was doing before which was
    # purely editing the node table itself. 
    # This plots the 

    # Get nodes with >0 edges and count how many occurences of each type there
    # are.
    result = node_table.sort(["edges"], ascending=False)
    connected_nodes = node_table[node_table["edges"] >0]
    connections = connected_nodes["edges"].tolist()
    connections.sort()
    counted = dict((i,connections.count(i)) for i in connections)    
    
    # convert dictionary to two lists for plotting
    x=[]
    y=[]
    for key, value in counted.iteritems():
        x.append(key)
        y.append(value)
    
    # Take log to make power law plot
    logx = [math.log(i,2) for i in x]
    logy = [math.log(j,2) for j in y]
    
    # Calculate trend line
    slope, intercept, r_value, p_value, std_err = linregress(logx,logy)
    ablineValues = []
    for i in logx:
      ablineValues.append(slope*i+intercept)
    
    # Plot this
    plot(logx,logy)
    plot(logx, ablineValues)
    show()
    return


def degree_abundance_presence_plots(node_list):
    # Plots node degree vs abundance and node degree vs presence.
    # Idea was to see if there is a relationship between these factors.
    # I can pretty comfortably say these factors are not important?
    
    # Next steps:
    # What about class or phylogenetic group? 
    # Maybe this could be cool. 

    connected = node_list[node_list["edges"]>0]
    connected["present"] = connected["av_abundance"]
    for taxa in connected.iterrows():
        connected["present"][connected["taxa"]==taxa[1][0]] =  np.count_nonzero(raw_data[raw_data["Taxon"] ==
        taxa[1][0]].iloc[0,1:-1])
        connected["av_abundance"][connected["taxa"] == taxa[1][0]] = float(raw_data[raw_data["Taxon"]==taxa[1][0]].mean(1))


    x = connected["edges"].tolist()
    y = connected["present"].tolist()
    z = connected["av_abundance"].tolist()
    #store_df(connected, path, "connected_nodes_df")

    logz = [math.log(i,10) for i in z]
    #scatter(x,logz)
    #show()
    return "these kinda suck"


#make_node_edge_tables(data, path)
#plot_abundance(node_list)
