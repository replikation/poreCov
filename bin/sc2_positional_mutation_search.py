################################################################################
## Module-import

import argparse
import pandas as pd
import numpy as np
import re


################################################################################
## function definition

def mutation_scan_setup(FASTA_FILE, MUTATION_LIST, SC2_GENES_START_END_NUCLEOTIDES_DF, FILE_MUTATIONS_DF):
  if FASTA_FILE != "Reference_genome.fasta":
    FILE_NAME = FASTA_FILE.replace('_filtered.consensus.fasta', '')
    FASTA_SEQUENCE = open(FASTA_FILE).readlines()[1]
    [mutation_scan(FILE_NAME, FASTA_SEQUENCE, SC2_GENES_START_END_NUCLEOTIDES_DF, FILE_MUTATIONS_DF, MUTATION) for MUTATION in MUTATION_LIST]
  
  return FILE_NAME, FILE_MUTATIONS_DF

def mutation_scan(FILE_NAME, FASTA_SEQUENCE, SC2_GENES_START_END_NUCLEOTIDES_DF, FILE_MUTATIONS_DF, MUTATION):
  if MUTATION not in list(FILE_MUTATIONS_DF.columns) or np.isnan(FILE_MUTATIONS_DF.loc[FILE_NAME, MUTATION]):
    MUTATION_DELIMITER_INDEX = MUTATION.find(":")
    if MUTATION_DELIMITER_INDEX != -1:
      MUTATION_GENE = MUTATION[0:MUTATION_DELIMITER_INDEX]
      MUTATION_CODON = re.sub("\D","", MUTATION[(MUTATION_DELIMITER_INDEX + 1):])
      if MUTATION_GENE in list(SC2_GENES_START_END_NUCLEOTIDES_DF['Gene']) and MUTATION_CODON != '':
        SIMILAR_MUTATIONS_LIST = [ELEMENT for ELEMENT in list(COMBINED_DF.columns) if MUTATION[0:(MUTATION_DELIMITER_INDEX + 1)] in ELEMENT if MUTATION_CODON == re.sub("\D","", ELEMENT[ELEMENT.index(':'):])]
        SIMILAR_MUTATIONS_VALUES_LIST = [COMBINED_DF.loc[FILE_NAME, SIMILAR_MUTATION] for SIMILAR_MUTATION in SIMILAR_MUTATIONS_LIST]
        if any(SIMILAR_MUTATIONS_VALUES_LIST): #check if any mutation already occured at the specified position
              FILE_MUTATIONS_DF.loc[FILE_NAME, MUTATION] = False
        
        else:
          MUTATION_GENE_START_NUCLEOTIDE = int(SC2_GENES_START_END_NUCLEOTIDES_DF[SC2_GENES_START_END_NUCLEOTIDES_DF["Gene"] == MUTATION_GENE]["Start_Nucleotide"])
          MUTATION_CODON_NR = int(MUTATION_CODON)
          CODON_START_NUCLEOTIDE = MUTATION_GENE_START_NUCLEOTIDE + ((MUTATION_CODON_NR - 1) * 3) #Position of the nucleotides of the mutated codon of each mutation
          CODON_NUCLEOTIDES_LIST = [FASTA_SEQUENCE[CODON_START_NUCLEOTIDE],FASTA_SEQUENCE[(CODON_START_NUCLEOTIDE + 1)],FASTA_SEQUENCE[(CODON_START_NUCLEOTIDE + 2)]]
          
          if "N" in CODON_NUCLEOTIDES_LIST:
            FILE_MUTATIONS_DF.loc[FILE_NAME, MUTATION] = 'Ambigous'
              
          else:
            FILE_MUTATIONS_DF.loc[FILE_NAME, MUTATION] = False
  
  return FILE_MUTATIONS_DF

def tsv_to_dataframe(TSV_FILE):
  if ".tsv" in TSV_FILE:
    OPENED_TSV_FILE = open(TSV_FILE).read()
    TSV_FILE_COLUMNS = OPENED_TSV_FILE.split('\n')[0].split("\t")
    TSV_FILE_VALUES = OPENED_TSV_FILE.split('\n')[1].split("\t")
    TSV_FILE_DF_COMPLETE = pd.DataFrame(TSV_FILE_VALUES, index = TSV_FILE_COLUMNS).transpose()
    COLUMNS = ['_id'] + TSV_FILE_DF_COMPLETE["aaSubstitutions"].iloc[0].split(',')
    VALUES = list([TSV_FILE_DF_COMPLETE["seqName"].iloc[0]]) + [True] * len(TSV_FILE_DF_COMPLETE["aaSubstitutions"].iloc[0].split(','))
    TSV_FILE_DF = pd.DataFrame(VALUES, index = COLUMNS).transpose().set_index('_id')

    return TSV_FILE_DF

################################################################################
## Initialization

parser = argparse.ArgumentParser(description = 'Scan SARS-CoV-2 fasta-file for sequence-state at mutation-sites from a list.')

parser.add_argument('-f', '--fasta', help = "Input fasta-file", required = True)
parser.add_argument('-m', '--mutationlist', help = "Input txt-file with SC2-mutations (one mutation/line)", required = True)
parser.add_argument('-n', '--name', help = "Name/ID of the sample", default = '')
parser.add_argument('-t', '--tsv', help = "Input tsv-file with clade & mutations", required = True)

#parsing:
arg = parser.parse_args()

#define arguments as variables:
FASTA_FILE = arg.fasta
MUTATION_LIST_FILE = arg.mutationlist
NAME = arg.name
TSV_FILE = arg.tsv


################################################################################
## Script

SC2_GENES_START_END_NUCLEOTIDES = {'Gene': ["E","M","N","ORF10","ORF1a","ORF1b","ORF3a","ORF6","ORF7a","ORF7b","ORF8","S","NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP11","NSP12","NSP13","NSP14","NSP15","NSP16"],
                                  'Start_Nucleotide': [26245,26523,28274,29558,266,13468,25393,27202,27394,27756,27894,21563,266,806,2720,8555,10055,10973,11843,12092,12686,13025,13442,13442,16237,18040,19621,20659],
                                  'End_Nucleotide': [26472,27191,29533,29674,13483,21555,26220,27387,27759,27887,28259,25384,805,2719,8554,10054,10972,11842,12091,12685,13024,13441,13480,16236,18039,19620,20658,21552]}
SC2_GENES_START_END_NUCLEOTIDES_DF = pd.DataFrame(data = SC2_GENES_START_END_NUCLEOTIDES)
#print(SC2_GENES_START_END_NUCLEOTIDES_DF)

MUTATION_LIST = open(MUTATION_LIST_FILE).read().split('\n')
#print(MUTATION_LIST)
tsv_to_dataframe(TSV_FILE)
#pd.set_option('display.expand_frame_repr', False)
#print(TSV_FILE_DF)

mutation_scan_setup(FASTA_FILE, MUTATION_LIST, SC2_GENES_START_END_NUCLEOTIDES_DF, TSV_FILE_DF)

if NAME != '':
  FILE_NAME = NAME

OUTPUT_FILE = FILE_NAME + "_positional_mutation_search.csv"
RESULT_CSV = TSV_FILE_DF.to_csv()
RESULT_FILE_CSV = open(OUTPUT_FILE, "w")
RESULT_FILE_CSV.write(result_csv)
RESULT_FILE_CSV.close()