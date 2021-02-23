#!/bin/bash
#Info: Creates report.csv for RKI from all pangolin.csv-files in the actual working dir.

INPUT_NAME=$1
OUTPUT_NAME=$2

echo "SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;PUBLICATION_STATUS;OWN_FASTA_ID" > $OUTPUT_NAME

while IFS=$'\t' read -r -a Array; do
    SENDING_LAB=$(echo ${SEQUENCING_LAB_ID})
    DATE_DRAW=$(echo "YYYYMMDD")
    SEQU_TYPE=$(echo "OXFORD_NANOPORE")
    SEQU_REASON=$(echo "X")
    SAMPLE_TYPE=$(echo "X")
    OWN_FASTA_ID="${Array[1]}"
    echo "${SENDING_LAB};${DATE_DRAW};${SEQU_TYPE};${SEQU_REASON};${SAMPLE_TYPE};N;${OWN_FASTA_ID}" >> $OUTPUT_NAME
done < ${INPUT_NAME}


