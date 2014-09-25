# Data formatting functions for generating and manipulating dataframes, getting OTU names from taxa ID strings etc.

def taxa_string_split(df, name_column, taxa_cat = [ "kingdom", "phylum", "class", 
                                                  "order", "family","genus"], depth=6):
    """ df[name_column] -> df[kingdom, phylum, class...] to desired depth
        This assumes you have a dataframe with one (or more) name columns. 
        By default it splits a taxa string in the QIIME format:
        "k__Bacteria; p__Verrucomicrobia; c__Opitutae; o__Opitutales; f__Opitutaceae;
        g__Opitutus; s__"
        into individual columns to make comparisons at many different phylogenetic levels possible.
         """
    for i in range(depth):
        splits = df[name_column].str.split().str.get(i)
        df[taxa_cat[i]] = splits

    return df


def dict_key_from_value(a_dict,a_value):
    """ returns first key that matches a dict value"""
    for key,value in a_dict.iteritems():
        if value == a_value:
            return key


def load_otu_sequences(rep_set, tax_assignments):
    """ Load files for sequences and names into a dataframe"""
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
