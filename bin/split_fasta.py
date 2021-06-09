#!/usr/bin/env python3
# SK

import os
import sys
import fnmatch

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)

def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

###

fasta_files = []
sequence_names = []
outfh = None

for fasta_file in sys.argv[1:]:

    if not os.path.exists(fasta_file):
        error(f'Not a file: {fasta_file}')

    log(f'Reading {fasta_file} ...')
    with open(fasta_file) as infh:
        for line in infh:
            if line.startswith('>'):

                # new sequence
                seq_name = line.strip().split()[0][1:]
                assert seq_name != '', f'Empty header in file: {fasta_file}'
                
                # sanitize
                seq_name = seq_name.replace('/', '_').replace(':', '_').replace('|','_')

                # handle duplicates
                if seq_name in sequence_names:
                    log(f'WARNING: Duplicate sequence name: {seq_name}')
                    # add number to seq_name, according to how often seq_name alrready appeared
                    duplicate_seq_name_pattern = f'{seq_name}_duplicate*'
                    seq_name += f'_duplicate{int(len([item for item in sequence_names if item == seq_name or fnmatch.fnmatch(item ,duplicate_seq_name_pattern)])):02d}'
                
                # save
                sequence_names.append(seq_name)

                # write to separate file
                if outfh is not None:
                    outfh.close()
                outfile = 'split_fasta/' + seq_name + '.fasta'
                log(f'Writing {outfile}')
                outfh = open(outfile, 'w')
                outfh.write(f'>{seq_name}\n')
                
            else:
                # write rest of lines (and fix windows line endings)
                outfh.write(line.replace('\r',''))

        # done for this file
        outfh.close()

log('Done.')
