#!/usr/bin/env bash

fasta_file=$1
txt_file=$2

grep ">" $fasta_file | grep ">" *.fasta | grep "complete genome" | cut -f2 -d "|" > .complete_genome_list.tmp

# create metadata file
printf "strain\tcountry\tdate\n" > clean_data.tsv
# create clean sequence file
rm "clean.references.fasta" 2>/dev/null  && touch "clean.references.fasta"

while IFS= read -r accession || [[ -n "$line" ]]; do

    entry=$(sed -n "/^ID   $accession;/,/\/\//p" $txt_file) 

    unset country date

    country=$(echo "$entry" | grep "^FT" | grep '/country="' | cut -f2 -d '"' |  tr ":|/ ,'" "_")
        if [ -z "$country" ]; then country="unknown"; fi
    date=$(echo "$entry" | grep "^FT" | grep ' /collection_date="' | cut -f2 -d '"')
        # skip entry if no date is available
        if [ -z "$date" ]; then continue; fi
    sequence=$(printf "$entry" | sed -n "/^SQ/,/\/\//p" | sed -e '2,$!d' -e '$d' | tr -d "0123456789 " | tr -d "\n" ) 

    # add data to metadatafile
    printf "${accession}-${country}-${date}\t${country}\t${date}\n" >> clean_data.tsv
    # exctract and add this sequence
    printf ">${accession}-${country}-${date}\n" >> clean.references.fasta
    printf "${sequence}\n" >> clean.references.fasta

done < .complete_genome_list.tmp