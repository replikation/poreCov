#!/usr/bin/env python3

################################################################################
## Module-import

import argparse
import pandas as pd
import numpy as np
import re


################################################################################
## function definition

def assign_to_dict(FILE_NAME, FILE_DICT):
  FILE_DICT[FILE_NAME] = pd.read_csv(FILE_NAME.strip("[],"))
  return FILE_DICT

################################################################################
## Initialization

parser = argparse.ArgumentParser(description = 'Combine output csv-files from "positional_mutation_scan.py" into one file.')

parser.add_argument('-f', '--files', nargs = '+', help = "Input csv-files", required = True)

#parsing:
arg = parser.parse_args()

#define arguments as variables:
FILES = arg.files


################################################################################
## Script

FILE_DICT = {}
[assign_to_dict(FILE_NAME, FILE_DICT) for FILE_NAME in FILES]

COMBINED_DF = pd.concat(FILE_DICT, axis = 0, ignore_index= True).set_index('_id')

RESULT_CSV = COMBINED_DF.to_csv()
RESULT_FILE_CSV = open("positional_mutation_summary.csv", "w")
RESULT_FILE_CSV.write(RESULT_CSV)
RESULT_FILE_CSV.close()