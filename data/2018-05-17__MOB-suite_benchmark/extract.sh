#!/bin/bash

# Extracts FASTA entry (header and sequence) from FASTA file, here MOB-suite_benchmark_seqs.fasta, 
# based on the specified accession number.



FASTA_FILE=/data/2018-05-17__MOB-suite_benchmark/MOB-suite_benchmark_seqs.fasta
ACC=$1
OUTPUT=$2


if [ -z "${ACC}" ] || [ -z "${OUTPUT}" ]; then

    echo -e "Error: Missing argument."
    echo -e "\nCorrect usage: extract.sh <accession number> <output file>"

else

    # determine line numbers (of the headers) of queried entry and the subsequent one
    LINE_NUMBERS=($(grep -n '^>' ${FASTA_FILE} | grep -E -A 1 ">${ACC}." | cut -d ':' -f1))
    
    # no entry with the specified accession number
    if [ "${#LINE_NUMBERS[@]}" == 0 ]; then
        
        echo "Error: There is no entry with accession ${ACC}. Please check ${FASTA_FILE}."    
    
    # accession number belongs to last entry, thus read to end of file
    elif [ "${#LINE_NUMBERS[@]}" == 1 ]; then
    
        echo -n "Extracting entry ${ACC} from ${FASTQ_FILE} to ${OUTPUT}..."    
        FIRST=${LINE_NUMBERS[0]}
        LAST=$(wc -l ${FASTA_FILE} | cut -d " " -f1)
        sed -n ${FIRST},${LAST}p ${FASTA_FILE} > ${OUTPUT}
        echo "DONE"
    
    # there is at least one entry after the queried one, 
    # thus read to the line before the header of that entry
    else
        
        echo -n "Extracting entry ${ACC} from ${FASTQ_FILE} to ${OUTPUT}..."        
        FIRST=${LINE_NUMBERS[0]}
        LAST=$((${LINE_NUMBERS[1]} - 1))
        sed -n ${FIRST},${LAST}p ${FASTA_FILE} > ${OUTPUT}
        echo "DONE"
        
    fi

fi