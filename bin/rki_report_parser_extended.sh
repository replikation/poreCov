#!/bin/bash
#Info: Creates report.csv for RKI from all pangolin.csv-files in the actual working dir.

INPUT_NAME=$1   # fastaname (array0), true value(array1)
OUTPUT_NAME=$2
EXTENDED_TABLE=$3

echo "SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;PUBLICATION_STATUS;OWN_FASTA_ID" > $OUTPUT_NAME

while IFS=$'\t' read -r -a Array; do
    SENDING_LAB=$(grep "${Array[0]}" ${EXTENDED_TABLE} | cut -f 2 -d ",")
    DATE_DRAW=$(grep "${Array[0]}" ${EXTENDED_TABLE} | cut -f 3 -d ",")
    SEQU_TECH=$(echo "OXFORD_NANOPORE")
    SEQU_REASON=$(grep "${Array[0]}" ${EXTENDED_TABLE} | cut -f 4 -d ",")
    SAMPLE_TYPE=$(grep "${Array[0]}" ${EXTENDED_TABLE} | cut -f 5 -d ",")
    OWN_FASTA_ID="${Array[0]}"
    echo "${SENDING_LAB};${DATE_DRAW};${SEQU_TECH};${SEQU_REASON};${SAMPLE_TYPE};N;${OWN_FASTA_ID}" >> $OUTPUT_NAME
done < ${INPUT_NAME}