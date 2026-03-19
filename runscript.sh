#!/bin/bash


echo "==============================================="
echo "Welcome to the gene finder!"
echo "==============================================="


echo "Choose an operation:"
echo "1. Find all genes from species 1 with at least one homolog in species 2."
echo "2. Find all human genes that have no homolog in mouse but at least one homolog in chimp."
echo "3. Find out if there are genes in orthologous groups that only contain mouse and rat proteins."
read -p "Enter option 1, 2 or 3: " option

if [ "$option" == "1" ]; then
    read -p "Enter species name 1: " species1
    read -p "Enter species name 2: " species2

    if [ "$species1" == "$species2" ]; then
        echo "Error: Species names must be different."
        exit 1
    fi
    echo "Running option 1..."
    homolog="homolog"
    python3 main.py "$species1" "$species2" "$homolog"

elif [ "$option" == "2" ]; then
    echo "Running option 2..."
    species3="Homo sapiens"
    species4="Mus musculus"
    species5="Pan troglodytes"
    python3 main.py "$species3" "$species4" "$species5"

elif [ "$option" == "3" ]; then
    echo "Running option 3..."
    species4="Mus musculus"
    species6="Rattus norvegicus"
    non="no homolog"
    python3 main.py "$species4" "$species6" "$non"

else
    echo "Invalid option. Please enter 1, 2, or 3."
    exit 1
fi