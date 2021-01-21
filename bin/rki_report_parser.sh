#!/bin/bash
#Info: Creates report.csv for RKI from all pangolin.csv-files in the actual working dir.

SENDING_LAB_ID=$1

echo "IMS_ID,SENDING_LAB,DATE_DRAW,SEQ_TYPE,SEQ_REASON,SAMPLE_TYPE,OWN_FASTA_ID" > rki_report.csv

for FILENAME in lineage*.csv; do
    IMS_ID=$(echo "IMS-00000-CVDP-00000")
    SENDING_LAB=$(echo "$SENDING_LAB_ID")       #In Process export "SENDING_LAB_ID" (given as Input/default when using [--rki]).
    DATE_DRAW=$(echo "YYYYMMDD")
    SEQU_TYPE=$(echo "OXFORD_NANOPORE")
    SEQU_REASON=$(echo "X")
    SAMPLE_TYPE=$(echo "X")
    OWN_FASTA_ID=$(tail -n+2 "$FILENAME" | rev | cut -f 6-10 -d "," | rev)
    echo $IMS_ID","$SENDING_LAB","$DATE_DRAW","$SEQU_TYPE","$SEQU_REASON","$SAMPLE_TYPE","$OWN_FASTA_ID >> rki_report.csv
done


