#!/bin/bash

# ==============================================================================
# Script Name: BLASTp_Db_Search.sh To search E.coli R388 T4SS Proteins Against pNOB8 Database
# Description: This script performs BLASTp searches for E.coli R388 T4SS 
#              proteins against a BLAST database created from the pNOB8 plasmid CDS sequences.
#              The script extracts CDS sequences from a GenBank file, creates a
#              BLAST database and executes BLASTp searches with specified proteins.
#
# References:
#   - BLAST® Command Line Applications User Manual [Internet]. Bethesda (MD): 
#     National Center for Biotechnology Information (US); 2008-. Available from:
#     https://www.ncbi.nlm.nih.gov/books/NBK279690/
#   - Websites used to help produce this script: https://bash.cyberciti.biz/guide/If..else..fi ; https://trstringer.com/python-in-shell-script/
#
# Requirements:
#   - Biopython installed for CDS extraction from GenBank files.
#   - BLASTtools command line tools installed for database creation and BLASTp searches.
# ==============================================================================


# Path to the GenBank files
genbank_path="/home/ib156/Desktop/Independant_research"

# Base path to the directories containing the protein FASTA files
protein_path="/home/ib156/Desktop/Independant_research/E.coli"

# Output directory for BLAST results
output_dir="/home/ib156/Desktop/Independant_research/results_w_qcov"
mkdir -p "$output_dir"

# Create an associative array for the protein files with their directories
declare -A proteins=(
    ["Relaxase_TrwC"]="Relaxase_TrwC/Q47673.fasta"
    ["VirB1_TrwN"]="VirB1_TrwN/Q6I6C7.fasta"
    ["VirB2_TrwL"]="VirB2_TrwL/O50328.fasta"
    ["VirB3_TrwM"]="VirB3_TrwM/O50329.fasta"
    ["VirB4_TrwK"]="VirB4_TrwK/O50330.fasta"
    ["VirB5_TrwJ"]="VirB5_TrwJ/O50331.fasta"
    ["VirB6_TrwI"]="VirB6_TrwI/O50333.fasta"
    ["VirB7_TrwH"]="VirB7_TrwH/O50334.fasta"
    ["VirB8_TrwG"]="VirB8_TrwG/O50335.fasta"
    ["VirB9_TrwF"]="VirB9_TrwF/O50336.fasta"
    ["VirB10_TrwE"]="VirB10_TrwE/O50337.fasta"
    ["VirB11_TrwD_ATPase"]="VirB11_TrwD_ATPase/O50338.fasta"
    ["VirD4_TrwB_ATPase"]="VirD4_TrwB_ATPase/Q04230.fasta"
)

# Path to the pNOB8 GenBank file, obtained from https://www.ncbi.nlm.nih.gov/nuccore/NC_006493
pNOB8_genbank_file="${genbank_path}/pNOB8.gb"

# Check if the pNOB8 GenBank file exists, to ensure directory is correct
if [ ! -f "$pNOB8_genbank_file" ]; then
    echo "GenBank file $pNOB8_genbank_file not found!"
    exit 1
fi

# Extract CDS sequences and create a FASTA file for pNOB8 containing these CDS sequences
python3 <<EOF    #Inspired from https://trstringer.com/python-in-shell-script/ for how to run python code within shell script
from Bio import SeqIO

pNOB8_fasta = 'pNOB8_cds.fasta'

with open(pNOB8_fasta, 'w') as pNOB8_handle:
    for record in SeqIO.parse('${pNOB8_genbank_file}', "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
                product = feature.qualifiers.get("product", ["unknown protein"])[0]
                translation = feature.qualifiers.get("translation", [""])[0]
                if translation:
                    fasta_entry = f">{protein_id} {product}\n{translation}\n"
                    pNOB8_handle.write(fasta_entry)
    print(f"CDS sequences have been extracted and written to {pNOB8_fasta}")
EOF

# Create BLAST database for pNOB8
if [ -f "pNOB8_cds.fasta" ] && [ -s "pNOB8_cds.fasta" ]; then
    makeblastdb -in pNOB8_cds.fasta -dbtype prot -out pNOB8_db
else
    echo "pNOB8_cds.fasta is missing or empty. Skipping pNOB8 database creation."
fi

# Perform BLAST search for each protein from E.coli's R388 T4SS against the pNOB8 database
for protein in "${!proteins[@]}"; do
    protein_file="${protein_path}/${proteins[$protein]}"
    if [ -f "$protein_file" ]; then
        if [ -f "pNOB8_db.pdb" ]; then
            blastp -query "$protein_file" -db pNOB8_db -out "${output_dir}/${protein}_vs_pNOB8.out" -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
        else
            echo "pNOB8_db.pdb not found! Skipping BLAST for $protein."
        fi
    else
        echo "File $protein_file not found!"
    fi    #source: https://bash.cyberciti.biz/guide/If..else..fi
done

echo "BLAST searches completed."

