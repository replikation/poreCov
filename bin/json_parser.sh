#!/usr/bin/env bash

SCRIPTNAME=$(basename -- "$0")

unset HASHID_INPUT PANGOLIN_INPUT NEXTCLADE_INPUT PRESIDENT_INPUT SEQU_LAB_ID PRIMER_INPUT
###############
##  Modules  ##
###############

usage()
  {
    echo "Parser script for PoreCov-workflow to generate MongoDb-report"
    echo " "
    echo "Usage:    $SCRIPTNAME [-i hashId ] [-p pangolin_result_file.csv] [-q president_result_file.tsv]"
    echo "          [-r sequencing_lab_id ] [-v primer-set] [-h help]"
    echo "Inputs:"
    echo -e "          [-i]    HashID of the sample"
    echo -e "          [-p]    Path to pangolin-result-file; .csv"
    echo -e "          [-n]    Path to nextclade-result-file; .tsv"
    echo -e "          [-q]    Path to president-result-file; .tsv"
    echo -e "          [-r]    ID of the sequencing-laboratory"
    echo -e "          [-v]    Used primer-set"
    exit;
  }

json_file_opening() {
    echo "{" > "$HASHID_INPUT"_mongodb_report.json
}

hashid_parsing() {
    echo "    \"_id\": {" >> "$HASHID_INPUT"_mongodb_report.json
    echo "        \"\$oid\": \"$HASHID_INPUT\"" >> "$HASHID_INPUT"_mongodb_report.json
    echo "    }," >> "$HASHID_INPUT"_mongodb_report.json
}

status_parsing() {
    echo "    \"Status\": \"analysed\"," >> "$HASHID_INPUT"_mongodb_report.json
}

lineage_parsing() {
    LINEAGE=$(tail -n +2 $PANGOLIN_INPUT | rev | cut -f 5 -d "," |rev)
    echo "    \"Lineage\": \"$LINEAGE\"," >> "$HASHID_INPUT"_mongodb_report.json
}

clade_parsing() {
    # clade is in row number 2 ($2 in the awk)
    CLADE=$(cat $NEXTCLADE_INPUT | awk -F "\t" '{print $2}' | grep -vw "clade")
    echo "    \"Clade\": \"$CLADE\"," >> "$HASHID_INPUT"_mongodb_report.json
}

mutation_parsing() {
    # AA mutations are in row number 17 ($17 in the awk)
    AA_MUTATIONS=$(cat $NEXTCLADE_INPUT | awk -F "\t" '{print $17}' | grep -vw "aaSubstitutions" |\
    tr "," "\n" | sed  's/$/": "true",/' | sed  's/^/        "/' | sed '$ s-.$--')

    echo "    \"Mutations\": {" >> "$HASHID_INPUT"_mongodb_report.json
    echo "$AA_MUTATIONS" >> "$HASHID_INPUT"_mongodb_report.json
    echo "    }," >> "$HASHID_INPUT"_mongodb_report.json
}

submitting_lab_parsing() {
    echo "    \"Submitting_Lab\": ," >> "$HASHID_INPUT"_mongodb_report.json
}

# if true remove line; else print number
sequ_lab_parsing() {
    echo "    \"Sequencing_Lab\": $SEQU_LAB_ID," >> "$HASHID_INPUT"_mongodb_report.json
}

primer_parsing() {
    echo "    \"Primer\": \"$PRIMER_INPUT\"," >> "$HASHID_INPUT"_mongodb_report.json
}

analysing_date_parsing() {
    ANALYSING_DATE=$(date -I | tr -d "-")
    echo "    \"Analysing_Date\": $ANALYSING_DATE," >> "$HASHID_INPUT"_mongodb_report.json
}

# if True write true else false
rki_valid_parsing() {
    RKI_VALID=$(tail -n +2 $PRESIDENT_INPUT | cut -f 16 -d $'\t' | tr '[:upper:]' '[:lower:]')
    echo "    \"RKI_Valid\": \"$RKI_VALID\"," >> "$HASHID_INPUT"_mongodb_report.json
}

rki_submit_parsing() {
    if [ $RKI_VALID == "false" ]; then
        echo "    \"RKI_Submit\": \"false\"," >> "$HASHID_INPUT"_mongodb_report.json
    fi
}

nucleotide_identity_parsing() {
    ACGT_NUCLEOTIDE_IDENTITY=$(tail -n +2 $PRESIDENT_INPUT | cut -f 1 -d $'\t')
    ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS=$(tail -n +2 $PRESIDENT_INPUT | cut -f 2 -d $'\t')
    if [ $RKI_VALID == "true" ]; then
        echo "    \"ACGT_Nucleotide_Identity\": $ACGT_NUCLEOTIDE_IDENTITY," >> "$HASHID_INPUT"_mongodb_report.json
        echo "    \"ACGT_Nucleotide_Identity_ignoring_Ns\": $ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS," >> "$HASHID_INPUT"_mongodb_report.json
    fi
}

ambiguous_bases_parsing() {
    AMBIGUOUS_BASES=$(tail -n +2 $PRESIDENT_INPUT | cut -f 8 -d $'\t')
    echo "    \"Ambiguous_Bases\": $AMBIGUOUS_BASES," >> "$HASHID_INPUT"_mongodb_report.json
}

query_length_parsing() {
    QUERY_LENGTH=$(tail -n +2 $PRESIDENT_INPUT | cut -f 13 -d $'\t')
    echo "    \"Query_Length\": $QUERY_LENGTH" >> "$HASHID_INPUT"_mongodb_report.json  #removed "," here, because itś the last entry
}

json_file_closing() {
    echo "}" >> "$HASHID_INPUT"_mongodb_report.json
}

#############################
###   Start of script    ####
#############################

while getopts 'i:p:n:q:r:v:h' flag; do
    case "${flag}" in
      i) HASHID_INPUT="${OPTARG}" ;;
      p) PANGOLIN_INPUT="${OPTARG}" ;;
      n) NEXTCLADE_INPUT="${OPTARG}" ;;
      q) PRESIDENT_INPUT="${OPTARG}" ;;
      r) SEQU_LAB_ID="${OPTARG}" ;;
      v) PRIMER_INPUT="${OPTARG}" ;;
      h) usage;;
      *) usage
         exit 1 ;;
    esac
done

# json head
json_file_opening; hashid_parsing
# inputfields
if [ ! -z "${SEQU_LAB_ID}" ]; then sequ_lab_parsing ; fi
if [ ! -z "${PRIMER_INPUT}" ]; then primer_parsing ; fi
status_parsing; lineage_parsing; clade_parsing; mutation_parsing;
analysing_date_parsing; rki_valid_parsing; rki_submit_parsing; nucleotide_identity_parsing;
ambiguous_bases_parsing; 
# last entry no comma
query_length_parsing
# json closure
json_file_closing