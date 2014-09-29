# associatio network functions

def testfun():
    print 'hello'
    return

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



