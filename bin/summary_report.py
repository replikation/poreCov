#!/usr/bin/env python3
'''
Generate a summary report for multiple samples
./summary_report.py
'''
# SK

import sys
import time
import argparse
import pandas as pd


def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\nAborting.\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

###

class SummaryReport():

    tabledata = None
    col_descriptions = []


    def validate_index(self, t_index):
        '''Assert that an index matches the tabledata index.'''
        if self.tabledata is None:
            error('Cannot validate_index when tabledata is None')
        if not sorted(self.tabledata.index) == sorted(t_index):
            error(f'Index mismatch: tabledata: {sorted(self.tabledata.index)} \n query: {sorted(t_index)}')


    def check_and_init_tabledata(self, t_index):
        '''If tabledata is None, initialize it now, then assert that an index matches the tabledata index.'''
        if self.tabledata is None:
            self.tabledata = pd.DataFrame(index=sorted(t_index))
            self.tabledata.columns.name = 'Sample'
        self.validate_index(t_index)


    def add_col_description(self, desc):
        self.col_descriptions.append(f'{desc}')


    def write_column_descriptions(self, filehandle):

        for desc in self.col_descriptions:
            filehandle.write(f'{desc}<br>\n')


    def write_html_report(self, filename='poreCov_summary_report.html'):
        '''Write the html report to a file'''

        htmlheader = '''<!DOCTYPE html><html><head>
        <title>poreCov Summary Report</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">

        <style>
        * {
            font-family:"Helvetica Neue",Helvetica,"Segoe UI",Arial,freesans,sans-serif
        }

        .content {
        max-width: 1600px;
        margin: auto;
        }

        table.tablestyle {
        background-color: #FFFFFF;
        width: 1600px;
        text-align: center;
        border-collapse: collapse;
        }
        table.tablestyle td, table.tablestyle th {
        border: 2px solid #8B8B8B;
        padding: 5px 5px;
        }
        table.tablestyle tbody td {
        font-size: 20px;
        color: #000000;
        }
        table.tablestyle tr:nth-child(even) {
        background: #E6F5FF;
        }
        table.tablestyle thead {
        background: #E6F5FF;
        }
        table.tablestyle thead th {
        font-size: 20px;
        font-weight: bold;
        color: #000000;
        text-align: center;
        }
        table.tablestyle tfoot td {
        font-size: 13px;
        }
        table.tablestyle tfoot .links {
        text-align: right;
        }
        table.tablestyle tfoot .links a{
        display: inline-block;
        background: #FFFFFF;
        color: #398AA4;
        padding: 2px 8px;
        border-radius: 5px;
        }
        </style>
        </head>

        <body>
        <div class="content">'''


        htmlfooter = '''
        </div>
        </body></html>
        '''

        htmltable = self.tabledata.to_html(classes=['tablestyle'])

        log(f'Writing report to {filename} ...')
        with open(filename, 'w') as outfh:
            outfh.write(htmlheader)

            outfh.write(f'''
            <h1 class="header" id="main-header">poreCov Summary Report</h1>
            Report created on: {time.asctime()}
            <h2 class="header" id="table-header">Sample overview</h2>
            ''')
            outfh.write(htmltable)

            self.write_column_descriptions(outfh)

            outfh.write(htmlfooter)



    ### functions to add columns

    def add_pangolin_results(self, pangolin_results):
        log(f'Adding Pangolin results ...')
        # column names:
        # taxon,lineage,probability,pangoLEARN_version,status,note
        res_data = pd.read_csv(pangolin_results, index_col=0)
        self.check_and_init_tabledata(res_data.index)

        res_data['lineage_prob'] = [f'{l} ({p})' for l,p in zip(res_data['lineage'], res_data['probability'])]
        self.tabledata['Lineage (probability)'] = res_data['lineage_prob']

        self.add_col_description('Lineage and probability were determined with Pangolin.')
            

    def add_president_results(self, president_results):
        log(f'Adding President results ...')
        # column names:
        # ACGT Nucleotide identity	ACGT Nucleotide identity (ignoring Ns)	ACGT Nucleotide identity (ignoring non-ACGTNs)
        #	Date	Matches	Mismatches	N_bases	Ngap	acgt_bases	file_in_query	file_in_refiupac_bases	length_query
        #	length_reference	non_upac_bases	qc_all_valid	qc_post_align_pass_threshold	qc_post_aligned	qc_post_aligned_all_valid
        # 	qc_valid_length	qc_valid_nucleotides	qc_valid_number_n	query_description	query_index	query_name	reference_name
        res_data = pd.read_csv(president_results, index_col=24, sep='\t')
        self.check_and_init_tabledata(res_data.index)

        self.tabledata['Nucleotide identity'] = res_data['ACGT Nucleotide identity']

        self.add_col_description('Nucleotide identity was determined with President.')

    
    def add_nextclade_results(self, nextclade_results):
        log(f'Adding Nextclade results ...')
        # column names:
        # seqName	clade	qc.overallScore	qc.overallStatus	totalGaps	totalInsertions	totalMissing	totalMutations
        # 	totalNonACGTNs	totalPcrPrimerChanges	substitutions	deletions	insertions	missing	nonACGTNs	pcrPrimerChanges
        # 	aaSubstitutions	totalAminoacidSubstitutions	aaDeletions	totalAminoacidDeletions	alignmentEnd	alignmentScore
        # 	alignmentStart	qc.missingData.missingDataThreshold	qc.missingData.score	qc.missingData.status
        # 	qc.missingData.totalMissing	qc.mixedSites.mixedSitesThreshold	qc.mixedSites.score	qc.mixedSites.status
        # 	qc.mixedSites.totalMixedSites	qc.privateMutations.cutoff	qc.privateMutations.excess	qc.privateMutations.score
        # 	qc.privateMutations.status	qc.privateMutations.total	qc.snpClusters.clusteredSNPs	qc.snpClust$rs.score
        # 	qc.snpClusters.status	qc.snpClusters.totalSNPs	errors
        res_data = pd.read_csv(nextclade_results, index_col=0, sep='\t')
        self.check_and_init_tabledata(res_data.index)

        res_data['mutations_formatted'] = [m.replace(',', ', ') for m in res_data['aaSubstitutions']]

        self.tabledata['Clade'] = res_data['clade']
        self.tabledata['Mutations'] = res_data['mutations_formatted']

        self.add_col_description('Clade and mutations were determined with Nextclade.')


    def add_kraken2_results(self, kraken2_results):
        log(f'Adding Kraken2 results ...')
        # column names:
        #sample,num_sarscov2,num_human
        res_data = pd.read_csv(kraken2_results, index_col=0)
        self.check_and_init_tabledata(res_data.index)

        res_data['total_reads'] = res_data['num_sarscov2'] + res_data['num_human']

        self.tabledata['%reads SARS-CoV-2'] = res_data['num_sarscov2'] / res_data['total_reads']
        self.tabledata['%reads human'] = res_data['num_human'] / res_data['total_reads']

        self.add_col_description('Read classification was determined with kraken2.')


###

if __name__ == '__main__':

    log('Started summary_report.py ...')

    parser = argparse.ArgumentParser(description='Generate a summary report for multiple samples poreCov')
    parser.add_argument("-p", "--pangolin_results", help="pangolin results")
    parser.add_argument("-n", "--nextclade_results", help="nextclade results")
    parser.add_argument("-q", "--president_results", help="president results")
    parser.add_argument("-k", "--kraken2_results", help="kraken2 results")
    args = parser.parse_args()


    # build report
    report = SummaryReport()

    # this determines the order of columns
    if args.kraken2_results:
        report.add_kraken2_results(args.kraken2_results)
    if args.pangolin_results:
        report.add_pangolin_results(args.pangolin_results)
    if args.president_results:
        report.add_president_results(args.president_results)
    if args.nextclade_results:
        report.add_nextclade_results(args.nextclade_results)

    
    report.write_html_report()


