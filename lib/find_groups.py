#!/usr/bin/env python3

import os
import re
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")

'''
The function find_oneTaxID() takes a species name and returns the species ID.

'''

def find_oneTaxID(species_name):
    species = []
    speciesID = None

    with open(os.path.join(DATA_DIR, 'eggnog4.species_list.txt'), 'r') as species_list:
        for line in species_list:
            if re.search(species_name, line):
                species = line.strip().split('\t')
    
    if species == []:
        print("Species name not found.")
        return None
    else:
        speciesID = species[1]
        return speciesID
    
'''
The function find_TaxIDs() takes two species names and returns the species IDs: 
[speciesID1, speciesID2]

'''

def find_TaxIDs(species_name1, species_name2):
    species_1 = []
    species_2 = []
    speciesID1 = None
    speciesID2 = None
    with open(os.path.join(DATA_DIR, 'eggnog4.species_list.txt'), 'r') as species_list:
        for line in species_list:
            if re.search(species_name1, line):
                species_1 = line.strip().split('\t')
            if re.search(species_name2, line):
                species_2 = line.strip().split('\t')
            if len(species_1) > 1 and len(species_2) > 1:
                break
            
    if species_1 == [] or species_2 == []:
        print("At least one species not found")
        return None
            
    speciesID1 = species_1[1]
    speciesID2 = species_2[1]
    
    return [speciesID1, speciesID2]

'''
The function findHomologProteins() takes two species IDs, optionally an input file containing meNOG members and a 
bolean value findHomolog. If findHomolog is true, the function creates prints the number of proteins in the first species 
have at least one homolog in the second species. 
If findHomolog is false, the function filters for all groups which only contain genes from species 1 and not species 2. (part used in 2a)
By default, the input file is: data/meNOG.members.tsv

'''

def findHomologProteins(speciesID1, speciesID2, speciesID = None, findHomolog = True):
    ID1 = str(speciesID1)
    ID2 = str(speciesID2)
    ID = str(speciesID)
    groups = []
    count = 0 
    with open(os.path.join(DATA_DIR, 'meNOG.members.tsv'), 'r') as mem_file:
        for line in mem_file:
            if findHomolog:
                if re.search(ID1, line) and re.search(ID2, line):
                    groups.append(line)
            else:
                if (re.search(ID1, line) 
                and not re.search(ID2, line)
                and re.search(ID, line)):
                    groups.append(line)
                    count += 1
                    
    if not groups and findHomolog:
        print(f"There's no gene in {ID1} that has at least one homolog in {ID2}.")
        return None
    if not groups and not findHomolog:
        print(f"There are no genes in {ID1} that have no homolog in {ID2}, but have at least one homolog in {ID}.")
    else:
        if findHomolog:
            print(f"{len(groups)} genes in {ID1} have at least one homolog in {ID2}.")
        else:
            print(f"{count} Homo sapiens genes do not have a homolog in mouse but have at least one in chimp.")  
            ids = []      
            with open("result_protein_ids.txt", "w") as f:
                f.write(f"The protein ids of the {count} Homo sapiens genes that do not have a homolog in mouse, but have homologs in chimp:")
                for line in groups:
                    last_col = line.strip().split('\t')[-1]
                    protein_id = last_col.split('.')[-1] # taking protein ids
                    ids.append(protein_id)
                f.write(f"{ids}\n")                     
    return groups

'''
The function details() takes as input what findHomologProteins() return, and that is the human genes in members.tsv that satisy the condition in 2a), and the meNOG.annotations.tsv file.
It selectes the corresponding functional categories, annotations and number of proteins for the human genes. In addition, it counts how many proteins are found in each functional category.
After it fetches all information required, it writes everything in a file called "results_annotations.tsv".
By default, the input file is: data/meNOG.annotations.tsv
'''

def details(groups, file = os.path.join(DATA_DIR, 'meNOG.annotations.tsv')):
    details_dict = {}
    with open(file, 'r') as details_file:
        for line in details_file:
            lines = line.strip().split('\t')
            group = lines[1]
            details_dict[group] = {
                "protein_count": lines[2],
                "functional_category": lines[4],
                "annotation": lines[5]
            }
            
    group_names = set()
    for line in groups:
        group_name = line.strip().split('\t')
        group_names.add(group_name[1])
        
    selected_dict = {}
    for name in group_names:
        if name in details_dict:
            selected_dict[name] = details_dict[name]
            
    functional_cat_counts = defaultdict(int)
    for info in selected_dict.values():
        functional_cat_counts[info["functional_category"]] += int(info["protein_count"])
        
    functional_cat_names = {}   
    with open(os.path.join(DATA_DIR, "eggnog4.functional_categories.txt"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') and ']' in line:
                letter = line[1]
                name = line.split(']', 1)[1].strip()
                functional_cat_names[letter] = name
    
    with open("results_annotations.tsv", "w") as out:
        out.write("GroupName\tProteinCount\tFunctionalCategory\tAnnotation\n")

        for name, info in selected_dict.items():
            out.write(
                f"{name}\t"
                f"{info['protein_count']}\t"
                f"{info['functional_category']}\t"
                f"{info['annotation']}\n"
            )

        out.write("\n# Functional category summary\n")
        out.write("Category\tProteinCount\n")

        for cat, count in sorted(functional_cat_counts.items()):
            cat_name = functional_cat_names.get(cat, "Multiple Functions")
            out.write(f"{cat} ({cat_name})\t{count}\n")
    
    return functional_cat_counts

'''
The function find_OrthoGroups() finds genes in orthologous groups that only contain rodent proteins (e.g. Mus Musculus, Rattus norvegicus).
It stores the orthologous group names, proteins ids and functional categories for the genes satisfying the condition in a file called "results_rodent.tsv". 
The input is the species IDs corresponding to Mus musculus and Rattus norvegicus, the meNOG.members.tsv file and the eggnog4.functional_categories.txt file.
'''
        
def find_OrthoGroups(species_ID1, species_ID2, file = os.path.join(DATA_DIR, 'meNOG.members.tsv'), file2 = os.path.join(DATA_DIR, 'eggnog4.functional_categories.txt')):
    ortho_dict = {}
    with open(file, 'r') as membersFile:
        for line in membersFile:
            cols = line.strip().split("\t")
            group_id = cols[1]
            functional_cat = cols[4]
            proteins_ = cols[5].split(",")
            species_ids_ = {p.split(".", 1)[0] for p in proteins_}
            
            if species_ids_ == {species_ID1, species_ID2}:
                protein_ids_ = [p.split(".", 1)[1] for p in proteins_]

                ortho_dict[group_id] = {
                    "protein_ids": protein_ids_,
                    "functional_category": functional_cat   
                }
                
    functional_cat_names = {}   
    with open(file2, 'r') as fe:
        for line in fe:
            line = line.strip()
            if line.startswith('[') and ']' in line:
                letter = line[1]
                name = line.split(']', 1)[1].strip()
                functional_cat_names[letter] = name
                
    with open("results_rodent.tsv", 'w') as out_file:
        out_file.write("GroupName\tProteinIDs\tFunctionalCategory\n")
        
        for name, info in ortho_dict.items():
            out_file.write(
                f"{name}\t"
                f"{info['protein_ids']}\t"
                f"{info['functional_category']}\n"
            )
        out_file.write("\n# Functional categories\n")
        out_file.write("Letter\tCategoryName\n")
        
        for name, info in ortho_dict.items():
            category_letter = info['functional_category'] 
            category_name = functional_cat_names.get(category_letter, "Combined name") 
            out_file.write(f"{category_letter} ({category_name})\n")
            
    return ortho_dict    
    
        
        
            
    