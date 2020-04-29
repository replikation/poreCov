#!/usr/bin/env bash

fasta_file=$1
maskBegin=$2 #'70'
maskEnd=$3 #'70'

pure_sequence=$(grep -v ">" $fasta_file | tr -d "\n" | tr -d " ") 



# check genome size
genome_size=$(echo $pure_sequence | wc -c )
# check not masked regions for amount of Ns
unspecific_bases=$(echo $pure_sequence | sed "s/^.\{$maskBegin\}//" | sed "s/.\{$maskEnd\}$//" | tr -d "ATGCatgc" | wc -c )


if (( $genome_size < 29000 ))
    then 
        echo "genom size for $fasta_file is to small; size: $genome_size bp" >> error_report_${fasta_file%.*}.txt
    else
        echo "all okay"
fi

if (( $unspecific_bases > 3 ));
    then 
        echo "To many unspecific bases for $fasta_file; Amount of N's = $unspecific_bases" >> error_report_${fasta_file%.*}.txt
    else
        echo "all okay"
fi