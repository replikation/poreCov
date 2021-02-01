#!/bin/bash
#Info: Creates report.csv for RKI from all pangolin.csv-files in the actual working dir.

SEQUENCING_LAB_ID=$1
INPUT_NAME=$2
OUTPUT_NAME=$3

echo "IMS_ID,SENDING_LAB,DATE_DRAW,SEQ_TYPE,SEQ_REASON,SAMPLE_TYPE,OWN_FASTA_ID" > $OUTPUT_NAME

while IFS=$'\t' read -r -a Array; do
    IMS_ID=$(echo "IMS-${SEQUENCING_LAB_ID}-CVDP-00000")
    SENDING_LAB=$(echo "00000")       #In Process export "SENDING_LAB_ID" (given as Input/default when using [--rki]).
    DATE_DRAW=$(echo "YYYYMMDD")
    SEQU_TYPE=$(echo "OXFORD_NANOPORE")
    SEQU_REASON=$(echo "X")
    SAMPLE_TYPE=$(echo "X")
    OWN_FASTA_ID="${Array[24]}"
    echo $IMS_ID","$SENDING_LAB","$DATE_DRAW","$SEQU_TYPE","$SEQU_REASON","$SAMPLE_TYPE","$OWN_FASTA_ID >> $OUTPUT_NAME
done < ${INPUT_NAME}


