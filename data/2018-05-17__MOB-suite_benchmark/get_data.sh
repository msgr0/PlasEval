#!/bin/bash

# Downloads the genome and plasmid sequences used in the evaluation of MOB-suite in
# "MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies"



LIST_FILE=SupplementalTable1.csv
FASTA_FILE=MOB-suite_benchmark_seqs.fasta


# read accessions into array
# assumptions: separator = semicolon, contains header row, accessions = first row
echo "Reading accessions from ${LIST_FILE}..."
ACCESSIONS=($(awk -F ";" 'NR > 1 {print $1}' ${LIST_FILE}))
echo "-> Number of accessions = ${#ACCESSIONS[@]}"

# download FASTA entries for all accessions and store them in one file
if [ -e ${FASTA_FILE} ]; then
    rm ${FASTA_FILE}
    echo "Deleted old version of ${FASTA_FILE}"
fi

echo "Downloading sequences..."
I=1
for ACC in ${ACCESSIONS[@]}; do
    echo -ne "\r${I} / ${#ACCESSIONS[@]}: ${ACC}"
    curl -s "https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" >> ${FASTA_FILE}
    I=$((I+1))
done

# summarise the download
N=$(grep '^>' ${FASTA_FILE} | wc -l)
L=$(grep -v '^>' ${FASTA_FILE} | tr -d '\n' | wc -c)
echo -e "\n-> Successfully downloaded ${N} sequences of total length ${L}"
echo "-> Output file: ${FASTA_FILE}"
