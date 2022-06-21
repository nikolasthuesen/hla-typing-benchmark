import pandas as pd
import re

#Create dict for conversion of deleted alleles (see deleted_alleles notebook)
deleted_df = pd.read_csv('../reference_data/Deleted_alleles.txt', comment='#')

def extract_allele_from_description(description):
    two_field_finder = re.search(r"(A|B|C|DRB1|DQB1)\*\d{2,3}:\d{2,4}", description)

    if two_field_finder != None:
        allele = two_field_finder.group(0)
    else:
        allele = ''

    return allele


#Only include alleles, which follow 2010 naming convention and have a proper new name (not e.g. just named by an error)
deleted_df = deleted_df[(deleted_df['Allele'].str.startswith(('A', 'B', 'C', 'DRB1', 'DQB1'))) & (~deleted_df['Allele'].str.startswith(('Cw')))]
deleted_df['allele_new'] = deleted_df.apply(lambda x: extract_allele_from_description(x['Description']), axis=1)
deleted_df = deleted_df[deleted_df['allele_new'] != '']

#Add two field naming of deleted allele:
deleted_df['allele_two_field'] = deleted_df.apply(lambda x: extract_allele_from_description(x['Allele']), axis=1)

#Only include alleles, that have a two field typing and where the two field typing differs from the new name
deleted_df = deleted_df[(deleted_df['allele_two_field'] != deleted_df['allele_new']) & (deleted_df['allele_two_field'] != '')]

#Ignore renaming of alleles in 3-field resolution, as the 1000G dataset and typing tools can't distinguish those
deleted_df = deleted_df[[len(split) < 3 for split in deleted_df['Allele'].str.split(':')]]

#Create dict for renaming
deleted_df = deleted_df.set_index('allele_two_field')['allele_new']

deleted_conversion_dict = deleted_df.to_dict()



#Function for converting an allele to one/two/three field resolution (disregarding any trailing letters - still unambiguous)
#For the genes A, B, C, DRB1 and DQB1 the first field always has two digits.

def convert_to_one_field(allele_high_res):
    import re
    one_field_finder = re.search(r"(A|B|C|DRB1|DQB1)\*\d{2}", allele_high_res)
    
    if one_field_finder != None:
        one_field_finder = one_field_finder.group(0)
    else:
        one_field_finder = None
    
    return one_field_finder    

    
def convert_to_two_field(allele_high_res):
    import re
    two_field_finder = re.search(r"(A|B|C|DRB1|DQB1)\*\d{2}:\d{2,4}", allele_high_res)
    
    if two_field_finder != None:
        allele_two_field = two_field_finder.group(0)
    else:
        allele_two_field = None
    
    #Finally, update the alleles, which have been changed/renamed            
    if allele_two_field in deleted_conversion_dict:
        allele_two_field = deleted_conversion_dict[allele_two_field]

    return allele_two_field    

def convert_to_three_field(allele_high_res):
    six_field_finder = re.search(r"(A|B|C|DRB1|DQB1)\*\d{2}:\d{2,4}:\d{2,4}", allele_high_res)
    
    if six_field_finder != None:
        allele_six_field = six_field_finder.group(0)
    else:
        allele_six_field = convert_to_two_field(allele_high_res)            
    
    return allele_six_field    


# Conversion to P-group and pseudosequence resolution

# Null allels are not listed in P group resolution .txt document. 
# However, these should still be mapped to their proper P group - if they are part of one. 
# The .txt file with G group resolution is therefore used to add these entries to p_group_dict
# Synonymous mutations are not grouped in G group, and these null alleles - even though they should be belong to a P group-
# are not grouped.

def make_p_group_dict(p_group_filepath='../reference_data/hla_nom_p.txt'):
    #Dict with the structure: {allele_in : allele_converted_to_p_group...}
    p_group_dict = dict()

    #Read the important results:
    with open(p_group_filepath, 'r') as infile:
        for line in infile:
            
            #If allele doesn't belong to a P group, add it to dict as key and value
            if ('/' not in line) and (line[0] != '#'):
                gene = line.split('*')[0]
                
                #Only register the valid alleles:
                if gene in ['A', 'B', 'C', 'DRB1', 'DQB1']:
                    p_group_entry = convert_to_two_field(gene + "*" + line.split(';')[-2])
                    
                    p_group_dict[p_group_entry] = p_group_entry
                    
            
            #If several alleles map to the same one, they are separated by a "/" - this indicates a P group
            if ('/' in line) and (line[0] != '#'):
                gene = line.split('*')[0]
                
                #Only register the valid alleles:
                if gene in ['A', 'B', 'C', 'DRB1', 'DQB1']:
                            
                    #Find the four field, P group resolution. The P group is found at the end of the line.
                    p_group_full = gene + "*" + line.split(';')[-1][:-1]
                    p_group_two_field = convert_to_two_field(p_group_full)
                
                    #Read the rest of the alleles and clean up the front and end part
                    synonymous_alleles = line.split('/')
                    synonymous_alleles[0] = synonymous_alleles[0].split(';')[1]
                    synonymous_alleles[-1] = synonymous_alleles[-1].split(';')[0]
                    

                    #Convert all alleles to four field resolution
                    for i in range(len(synonymous_alleles)):
                        synonymous_alleles[i] = gene + "*" + synonymous_alleles[i]
                        synonymous_alleles[i] = convert_to_two_field(synonymous_alleles[i])

                    #Remove duplicates when converting to four field:
                    synonymous_allels_unique_two_field = list(set(synonymous_alleles))
                    
                    #Add key in dict for each of the unique entries:
                    for synonymous_allele in synonymous_allels_unique_two_field:
                                                
                        p_group_dict[synonymous_allele] = convert_to_two_field(p_group_full)

    return p_group_dict


#Make P type conversion function using p_group_dict
def convert_to_p_group(allele, p_group_dict=''):
    if p_group_dict == '':
        p_group_dict = make_p_group_dict()

    #Start by converting to two field:
    allele_two_field = convert_to_two_field(allele)
    
    #Find corresponding P-type if it exists. 
    if allele_two_field in p_group_dict:
        allele_two_field_p_group = p_group_dict[allele_two_field]
    
    #Failsafe: if allele isn't found in p_group_dict (which it should be), then just use two_field_resolution
    else:
        allele_two_field_p_group = allele_two_field
        #print("Allele was not found in p_group_dict", allele_two_field)
       
        
    return allele_two_field_p_group



#Make dict for evaxion-group conversion
def make_e_group_dict(e_group_filepath='../reference_data/classic.mhc_seqs.tsv'):
    hla_loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1']

    mhc_pseudo_df = pd.read_csv(e_group_filepath, sep='\t')[['allele', 'mhc_seq']]
    mhc_pseudo_df = mhc_pseudo_df[mhc_pseudo_df['allele'].str.startswith(tuple(hla_loci))]
    mhc_pseudo_df['allele'] = mhc_pseudo_df['allele'].str.replace('HLA-', '')
    mhc_pseudo_df = mhc_pseudo_df.set_index('allele')
    e_group_dict = mhc_pseudo_df['mhc_seq'].to_dict()
    
    return e_group_dict
            


#Function for converting to e-group resolution:
def convert_to_e_group(allele, e_group_dict='', p_group_dict=''):
    if e_group_dict == '':
        e_group_dict = make_e_group_dict()
    if p_group_dict == '':
        p_group_dict = make_p_group_dict()

    #Start by converting to P group:
    allele_p_group = convert_to_p_group(allele, p_group_dict)
    
    #Find corresponding e-type if it exists
    if allele_p_group in e_group_dict:
        allele_e_group = e_group_dict[allele_p_group]
    
    #Else use P group
    else:
        # print("Using P group as backup")
        # print(allele)
        allele_e_group = allele_p_group
        
    return allele_e_group

