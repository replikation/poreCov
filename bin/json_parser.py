#!/usr/bin/env python3

################################################################################
## Module-import

import pandas as pd
import numpy as np
import os
import sys


################################################################################
## Initialization
import argparse
parser = argparse.ArgumentParser(description = 'Create json-file for upload to MongoDB from different result-files.')

parser.add_argument('-i', '--id', help = "Input hash-Id", required = True)
parser.add_argument('-n', '--nextclade', help = "Input Nextclade-file", required = True)
parser.add_argument('-o', '--output', help = "Output-directory", default = os.getcwd())
parser.add_argument('-p', '--pangolin', help = "Input Pangolin-file", required = True)
parser.add_argument('-q', '--president', help = "Input President-file", required = True)
parser.add_argument('-s', '--seqlab', help = "ID of the sequencing-laboratory", default = 'no ID given')
parser.add_argument('-t', '--sublab', help = "ID of the submitting-laboratory", default = 'no ID given')
parser.add_argument('-v', '--primer', help = "Used primer-set", default = '')

#parsing:
arg = parser.parse_args()

#define arguments as variables:
HASHID_INPUT = arg.id
NEXTCLADE_INPUT = arg.nextclade
OUTPUT_DIR = arg.output
PANGOLIN_INPUT = arg.pangolin
PRESIDENT_INPUT = arg.president
PRIMER_INPUT = arg.primer
SEQLAB_ID = arg.seqlab
SUBLAB_ID = arg.sublab


################################################################################
## Set up results-directory & change working-directory if necessary

#check if output-directory exists & create a Results-directory in it:
OUTPUT_PATH = str(OUTPUT_DIR) + '/'
os.makedirs(OUTPUT_PATH, exist_ok = True)

#check if output-flag was used:
if OUTPUT_PATH != os.getcwd():
        
    #change working directory to the path
    os.chdir(OUTPUT_PATH)


################################################################################
## Dataframe-creation

DF_NEXTCLADE = pd.read_csv(NEXTCLADE_INPUT, sep = '\t')
DF_PANGOLIN = pd.read_csv(PANGOLIN_INPUT)
DF_PRESIDENT = pd.read_csv(PRESIDENT_INPUT, sep = '\t')


################################################################################
## Parsing-functions

OUTPUT_FILE_NAME = str(HASHID_INPUT) + "_mongodb_report.json"

def json_file_opening(OUTPUT_FILE_NAME):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "w")
	RESULT_FILE.write("{\n")
	RESULT_FILE.close()
	return RESULT_FILE

def hashid_parsing(OUTPUT_FILE_NAME, HASHID_INPUT):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"_id\": {\n")
	RESULT_FILE.write("        \"$oid\": \"%s\"\n" %(HASHID_INPUT))
	RESULT_FILE.write("    },\n")
	RESULT_FILE.close()
	return RESULT_FILE

def status_parsing(OUTPUT_FILE_NAME):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Status\": \"analysed\",\n")
	RESULT_FILE.close()
	return RESULT_FILE

def lineage_parsing(OUTPUT_FILE_NAME, DF_PANGOLIN):
	LINEAGE = DF_PANGOLIN['lineage'].values[0]
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Lineage\": \"%s\",\n" %(LINEAGE))
	RESULT_FILE.close()
	return RESULT_FILE

def clade_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE):
	CLADE = DF_NEXTCLADE['clade'].values[0]
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Clade\": \"%s\",\n" %(CLADE))
	RESULT_FILE.close()
	return RESULT_FILE

def mutation_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE):
	AA_MUTATIONS_LIST = str(DF_NEXTCLADE['aaSubstitutions'].values[0]).split(',')
	if AA_MUTATIONS_LIST[0] == 'nan':
		AA_MUTATIONS_LIST[0] = 'no_mutations'
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Mutations\": {\n")
	[RESULT_FILE.write("        \"%s\": \"true\",\n" %(AA_MUTATION)) if AA_MUTATION != AA_MUTATIONS_LIST[-1] else RESULT_FILE.write("        \"%s\": \"true\"\n" %(AA_MUTATION)) for AA_MUTATION in AA_MUTATIONS_LIST]
	RESULT_FILE.write("    },\n")
	RESULT_FILE.close()
	return RESULT_FILE

def deletion_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE):
	AA_DELETIONS_LIST = str(DF_NEXTCLADE['aaDeletions'].values[0]).split(',')
	if AA_DELETIONS_LIST[0] == 'nan':
		AA_DELETIONS_LIST[0] = 'no_mutations'
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Deletions\": {\n")
	[RESULT_FILE.write("        \"%s\": \"true\",\n" %(AA_DELETION)) if AA_DELETION != AA_DELETIONS_LIST[-1] else RESULT_FILE.write("        \"%s\": \"true\"\n" %(AA_DELETION)) for AA_DELETION in AA_DELETIONS_LIST]
	RESULT_FILE.write("    },\n")
	RESULT_FILE.close()
	return RESULT_FILE

def submitting_lab_parsing(OUTPUT_FILE_NAME, SUBLAB_ID):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Submitting_Lab\": ,\n")
	RESULT_FILE.close()
	return RESULT_FILE

# if true remove line; else print number
def sequ_lab_parsing(OUTPUT_FILE_NAME, SEQLAB_ID):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Sequencing_Lab\": %s,\n" %(SEQLAB_ID))
	RESULT_FILE.close()
	return RESULT_FILE

def primer_parsing(OUTPUT_FILE_NAME, PRIMER_INPUT):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Primer\": \"%s\",\n" %(PRIMER_INPUT))
	RESULT_FILE.close()
	return RESULT_FILE

def analysing_date_parsing(OUTPUT_FILE_NAME):
	DATE = os.popen('date -I | tr -d "-" |tr -d "\n"')
	ANALYSING_DATE = DATE.read()
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Analysing_Date\": %s,\n" %(ANALYSING_DATE))
	RESULT_FILE.close()
	return RESULT_FILE

# if True write true else false
def rki_valid_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT):
	RKI_VALID = str(DF_PRESIDENT['qc_all_valid'].values[0]).lower()
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"RKI_Valid\": \"%s\",\n" %(RKI_VALID))
	RESULT_FILE.close()
	return RESULT_FILE, RKI_VALID

def rki_submit_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT):
	RKI_VALID = str(DF_PRESIDENT['qc_all_valid'].values[0]).lower()
	if RKI_VALID == "false":
		RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
		RESULT_FILE.write("    \"RKI_Submit\": \"false\",\n")
		RESULT_FILE.close()
		return RESULT_FILE

def nucleotide_identity_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT):
	RKI_VALID = str(DF_PRESIDENT['qc_all_valid'].values[0]).lower()
	if RKI_VALID == "true":
		ACGT_NUCLEOTIDE_IDENTITY = DF_PRESIDENT['ACGT Nucleotide identity'].values[0]
		ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS = DF_PRESIDENT['ACGT Nucleotide identity (ignoring Ns)'].values[0]
		RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
		RESULT_FILE.write("    \"ACGT_Nucleotide_Identity\": %s,\n" %(ACGT_NUCLEOTIDE_IDENTITY))
		RESULT_FILE.write("    \"ACGT_Nucleotide_Identity_ignoring_Ns\": %s,\n" %(ACGT_NUCLEOTIDE_IDENTITY_IGNORING_NS))
		RESULT_FILE.close()
		return RESULT_FILE

def ambiguous_bases_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT):
	AMBIGUOUS_BASES = DF_PRESIDENT['N_bases'].values[0]
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Ambiguous_Bases\": %s,\n" %(AMBIGUOUS_BASES))
	RESULT_FILE.close()
	return RESULT_FILE

def query_length_parsing(OUTPUT_FILE_NAME):
	QUERY_LENGTH = DF_PRESIDENT['length_query'].values[0]
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("    \"Query_Length\": %s\n" %(QUERY_LENGTH)) #removed "," here, because itś the last entry
	RESULT_FILE.close()
	return RESULT_FILE

def json_file_closing(OUTPUT_FILE_NAME):
	RESULT_FILE = open(OUTPUT_FILE_NAME, "a")
	RESULT_FILE.write("}")
	RESULT_FILE.close()
	return RESULT_FILE


################################################################################
## Function calls

json_file_opening(OUTPUT_FILE_NAME)
hashid_parsing(OUTPUT_FILE_NAME, HASHID_INPUT)
status_parsing(OUTPUT_FILE_NAME)
lineage_parsing(OUTPUT_FILE_NAME, DF_PANGOLIN)
clade_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE)
mutation_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE)
deletion_parsing(OUTPUT_FILE_NAME, DF_NEXTCLADE)
#sequ_lab_parsing(OUTPUT_FILE_NAME, SEQLAB_ID)
#submitting_lab_parsing(OUTPUT_FILE_NAME, SUBLAB_ID)
if PRIMER_INPUT != "":
	primer_parsing(OUTPUT_FILE_NAME, PRIMER_INPUT)
analysing_date_parsing(OUTPUT_FILE_NAME)
rki_valid_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT)
rki_submit_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT)
nucleotide_identity_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT)
ambiguous_bases_parsing(OUTPUT_FILE_NAME, DF_PRESIDENT)
query_length_parsing(OUTPUT_FILE_NAME)
json_file_closing(OUTPUT_FILE_NAME)