#!/usr/bin/env bash

name=$1
file_input=$2

mkdir -p ${name}_contigs/

while read line
   do
      if [[ ${line:0:1} == '>' ]]
      then
        outfile=$(printf "${line#>}" | tr "|" "_" | tr -d "\r" | tr "/" "_" | tr " " "_" )
        echo "${line}" | tr -d "\r" > ${name}_contigs/${outfile}.fa
      else
        echo "${line}" | tr -d "\r" >> ${name}_contigs/${outfile}.fa
      fi
   done < ${file_input}