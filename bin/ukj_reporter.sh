#!/usr/bin/env bash

SCRIPTNAME=$(basename -- "$0")

###############
##  Modules  ##
###############

usage()
  {
    echo "Parser script for PoreCov-workflow to generate MongoDb-report"
    echo " "
    echo "Usage:    $SCRIPTNAME [-i hashId ] [-p pangolin_result_file.csv] [-q president_result_file.tsv]"
    echo "          [-r seqeuncing_lab_id ] [-v primer-set] [-h help]"
    echo "Inputs:"
    echo -e "          [-i]    HashID of the sample"
    echo -e "          [-p]    Path to pangolin-result-file; .csv"
    echo -e "          [-q]    Path to president-result-file; .tsv"
    echo -e "          [-r]    ID of the sequencing-laboratory"
    echo -e "          [-v]    Used primer-set"
    exit;
  }

json_file_opening() {
    echo "{" > "$HASHID_INPUT"_mongodb_report.json
}

hashid_parsing() {
    echo '    "_id": {' >> "$HASHID_INPUT"_mongodb_report.json
    echo '        "$oid": "'"$HASHID_INPUT"'"' >> "$HASHID_INPUT"_mongodb_report.json
    echo '    },' >> "$HASHID_INPUT"_mongodb_report.json
}

status_parsing() {
    echo '    "Status": "analysed",' >> "$HASHID_INPUT"_mongodb_report.json
}

lineage_parsing() {
    LINEAGE=$(tail -n +2 $PANGOLIN_INPUT | rev | cut -f 5 -d "," |rev)
    echo '    "Lineage": "'"$LINEAGE"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

submitting_lab_parsing() {
    echo '    "Submitting_Lab": "",' >> "$HASHID_INPUT"_mongodb_report.json
}

sequ_lab_parsing() {
    echo '    "Sequencing_Lab": "'"$SEQU_LAB_ID"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

primer_parsing() {
    echo '    "Primer": "'"$PRIMER_INPUT"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

analysing_date_parsing() {
    ANALYSING_DATE=$(tail -n +2 $PRESIDENT_INPUT | cut -f 13 -d $'\t')
    echo '    "Analysing_Date": "'"$ANALYSING_DATE"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

rki_valid_parsing() {
    RKI_VALID=$(tail -n +2 $PRESIDENT_INPUT | cut -f 2 -d $'\t')>> "$HASHID_INPUT"_mongodb_report.json
    echo '    "RKI_Valid": "'"$RKI_VALID"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

rki_submit_parsing() {
    if [ $RKI_VALID == "False" ]; then
        echo '    "RKI_Submit": "False",' >> "$HASHID_INPUT"_mongodb_report.json
    fi
}

nucleotide_identity_parsing() {
    ACGT_NUCLEOTIDE_IDENTITY=$(tail -n +2 $PRESIDENT_INPUT | cut -f 3 -d $'\t')
    ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS=$(tail -n +2 $PRESIDENT_INPUT | cut -f 4 -d $'\t')
    echo '    "ACGT_Nucleotide_Identity": "'"$ACGT_NUCLEOTIDE_IDENTITY"'",' >> "$HASHID_INPUT"_mongodb_report.json
    echo '    "ACGT_Nucleotide_Identity_ignoring_Ns": "'"$ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

ambiguous_bases_parsing() {
    AMBIGUOUS_BASES=$(tail -n +2 $PRESIDENT_INPUT | cut -f 6 -d $'\t')
    echo '    "Ambiguous_Bases": "'"$AMBIGUOUS_BASES"'",' >> "$HASHID_INPUT"_mongodb_report.json
}

query_length_parsing() {
    QUERY_LENGTH=$(tail -n +2 $PRESIDENT_INPUT | cut -f 7 -d $'\t')
    echo '    "Query_Length": "'"$QUERY_LENGTH"'"' >> "$HASHID_INPUT"_mongodb_report.json  #removed "," here, because itś the last entry
}

json_file_closing() {
    echo "}" >> "$HASHID_INPUT"_mongodb_report.json
}

#############################
###   Start of script    ####
#############################

while getopts 'i:p:q:r:v:h' flag; do
    case "${flag}" in
      i) HASHID_INPUT="${OPTARG}" ;;
      p) PANGOLIN_INPUT="${OPTARG}" ;;
      q) PRESIDENT_INPUT="${OPTARG}" ;;
      r) SEQU_LAB_ID="${OPTARG}" ;;
      v) PRIMER_INPUT="${OPTARG}" ;;
      h) usage;;
      *) usage
         exit 1 ;;
    esac
done

# json head
json_file_opening
# inputfields
hashid_parsing; status_parsing; lineage_parsing; submitting_lab_parsing; sequ_lab_parsing
primer_parsing; analysing_date_parsing; rki_valid_parsing; rki_submit_parsing; nucleotide_identity_parsing
ambiguous_bases_parsing; query_length_parsing
# json closure
json_file_closing