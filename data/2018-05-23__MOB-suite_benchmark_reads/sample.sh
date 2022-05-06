#!/bin/bash

# Provides easy access to a sample by providing information on it or its FASTQ data.
# It also allows to extract the accession numbers of chromosomes and plasmids of a sample.
# Samples are accessed via their identifier (first column in samples.csv).
#
# Usage:
# see help message below


SAMPLES_FILE=/data/2018-05-23__MOB-suite_benchmark_reads/samples.csv
FASTQ_DIR=/data/2018-05-23__MOB-suite_benchmark_reads/fastq-split
CMD=$1
SID=$2
OUTDIR=$3

SRA=$2 # alternative id variable when an SRA accession is used as input (just for better readability)


# show details of the sample
if [ ${CMD} == "view" ] && [ ! -z "${SID}" ]; then

    SAMPLE=$(grep -E "^${SID};" ${SAMPLES_FILE})

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample with identifier ${SID}. Please check ${SAMPLES_FILE}."
    else
        echo ${SAMPLE} | awk -F ';' '{print "Sample ID:\t" $1; print "Organism name:\t" $2; print "Sample name:\t" $3; print "SRA accession:\t" $4; print "Chromosome:\t" $5; print "Plasmids:\t" $6}'
    fi

# show SRA accession
elif [ ${CMD} == "sra" ] && [ ! -z "${SID}" ]; then

    SAMPLE=$(grep -E "^${SID};" ${SAMPLES_FILE})

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample with identifier ${SID}. Please check ${SAMPLES_FILE}."
    else
        echo ${SAMPLE} | awk -F ';' '{print $4}'
    fi

# show chromosome accession
elif [ ${CMD} == "chr" ] && [ ! -z "${SID}" ]; then

    SAMPLE=$(grep -E "^${SID};" ${SAMPLES_FILE})

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample with identifier ${SID}. Please check ${SAMPLES_FILE}."
    else
        echo ${SAMPLE} | awk -F ';' '{print $5}'
    fi

# show plasmid accession(s)
elif [ ${CMD} == "pla" ] && [ ! -z "${SID}" ]; then

    SAMPLE=$(grep -E "^${SID};" ${SAMPLES_FILE})

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample with identifier ${SID}. Please check ${SAMPLES_FILE}."
    else
        echo ${SAMPLE} | awk -F ';' '{print $6}'
    fi

# extract FASTQ file of the sample
elif [ ${CMD} == "extract" ] && [ ! -z "${SID}" ] && [ ! -z "${OUTDIR}" ]; then

    OUTDIR=$(readlink -f $OUTDIR)

    SAMPLE=$(grep -E "^${SID};" ${SAMPLES_FILE})
    SRA=$(echo "${SAMPLE}" | cut -d ";" -f4)

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample with identifier ${SID}. Please check ${SAMPLES_FILE}."
    else

        FILE=${FASTQ_DIR}/${SRA}.fastq.bz2
        if [ -e ${FILE} ]; then
            echo -n "Extracting sample ${SID} from ${FILE} to ${OUTDIR}..."
            bzip2 -dkc "${FILE}" > "${OUTDIR}/${SRA}.fastq"
            echo "DONE"
        fi

        FILE=${FASTQ_DIR}/${SRA}_1.fastq.bz2
        if [ -e ${FILE} ]; then
            echo -n "Extracting sample ${SID} from ${FILE} to ${OUTDIR}..."
            bzip2 -dkc "${FILE}" > "${OUTDIR}/${SRA}_1.fastq"
            echo "DONE"
        fi

        FILE=${FASTQ_DIR}/${SRA}_2.fastq.bz2
        if [ -e ${FILE} ]; then
            echo -n "Extracting sample ${SID} from ${FILE} to ${OUTDIR}..."
            bzip2 -dkc "${FILE}" > "${OUTDIR}/${SRA}_2.fastq"
            echo "DONE"
        fi

    fi

# show sample id
elif [ ${CMD} == "sid" ] && [ ! -z "${SRA}" ]; then

    SAMPLE=$(grep -E ";{SRA};" ${SAMPLES_FILE})

    if [ -z "${SAMPLE}" ]; then
        echo "Error: There is no sample associated with SRA accession ${SRA}. Please check ${SAMPLES_FILE}."
    else
        echo ${SAMPLE} | cut -d ';' -f1
    fi

# print help
else

    echo -e "Error: No / unknown command or missing parameter(s).\n\nPossible commands:"
    echo -e "(1) sample.sh view <sample id>\n\tShows information associated with the sample identifier."
    echo -e "(2) sample.sh extract <sample id> <output directory>\n\tExtracts the sample into the specified directory."
    echo -e "(3) sample.sh sra <sample id>\n\tDetermines the SRA accession number of the sample."
    echo -e "(4) sample.sh chr <sample id>\n\tDetermines the accession numbers of the chromosome(s) in the sample."
    echo -e "(5) sample.sh pla <sample id>\n\tDetermines the accession numbers of the plasmid(s) in the sample."
    echo -e "(6) sample.sh sid <sra accession>\n\tDetermines the sample id associated with the SRA accession number."

fi

