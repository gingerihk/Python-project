#!/usr/bin/env python3
import sys
import os
import lib.find_groups as fg

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

arguments = sys.argv[1:]
if len(arguments) != 3:
    print(f"Error: main.py requires 3 arguments, {len(arguments)} are given.") 
    sys.exit()
    
species_name1 = arguments[0]
species_name2 = arguments[1]
if arguments[2] == 'homolog':
    homolog = True
if arguments[2] == 'Pan troglodytes':
    species_name = arguments[2]
    homolog = False
if arguments[2] == 'no homolog':
    homolog = False

speciesID1, speciesID2 = fg.find_TaxIDs(species_name1, species_name2)
if speciesID1 == None or speciesID2 == None:
    print("Try again with different species name.")
    sys.exit()

print(f"The Species ID for {species_name1} is: {speciesID1}.")
print(f"The Species ID for {species_name2} is: {speciesID2}.")

# Task 1: program execution: python3 main.py <species_name1> <species_name2> 'homolog'
if homolog:
    print(f"Searching for genes from {species_name1} which have at least one homolog in {species_name2}.")
    groups = fg.findHomologProteins(speciesID1, speciesID2)
    if groups == None:
        sys.exit()
        
# Task 2: program execution: python3 main.py 'Homo sapiens' 'Mus musculus' 'Pan troglodytes'
if not homolog:
    if arguments[0] == 'Homo sapiens':
        speciesID = fg.find_oneTaxID(species_name)
        print(f"The Species ID for {species_name} is: {speciesID}.")
        print(f"Searching for human genes that have no homolog in mouse, but at least one homolog in chimp.")
        speciesID1, speciesID2 = fg.find_TaxIDs(species_name1, species_name2)
        print(f"Writing their protein IDs in 'result_protein_ids.txt'.")
        outputFile = os.path.join(RESULTS_DIR, 'result_protein_ids.txt')
        if os.path.exists(outputFile):
            print(f"SKIP. {outputFile} already exists.")
        groups = fg.findHomologProteins(speciesID1, speciesID2, speciesID, findHomolog = False)
        print(f"Searching for and writing details (annotations, functional category, etc) corresponding to the genes of interest in 'results_annotations.tsv'.")
        outputFile_1 = os.path.join(RESULTS_DIR, 'results_annotations.tsv')
        if os.path.exists(outputFile_1):
            print(f"SKIP. {outputFile_1} already exists.")
        details_ = fg.details(groups)
        if groups == None:
            sys.exit()
            
# Task 3: program execution: python3 main.py 'Mus musculus' 'Rattus norvegicus' 'no homolog'  
    if arguments[0] == 'Mus musculus':
            species_ID_1 = fg.find_oneTaxID(species_name1)
            species_ID_2 = fg.find_oneTaxID(species_name2)
            if species_ID_1 == None or species_ID_2 == None:
                print("Try again with 'Mus musculus' and 'Rattus norvegicus'.")
                sys.exit()
            
            outputFile_2 = os.path.join(RESULTS_DIR,'results_rodent.tsv') 
            if os.path.exists(outputFile_2):
                print(f"Genes already found! SKIP. {outputFile_2} already exists.")
            else:
                print(f"Checking if there are genes in orthologous groups that only contain mouse and rat proteins")
                ortho_groups_ = fg.find_OrthoGroups(species_ID_1, species_ID_2)
                if ortho_groups_ == None:
                        sys.exit()
                print(f"Genes found! Results written in 'results_rodent.tsv'.")
            
            
                