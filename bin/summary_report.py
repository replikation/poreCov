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

    report_time = time.localtime()
    report_name = f'poreCov_summary_report_{time.strftime("%Y-%m-%d--%H-%M-%S", report_time)}'
    output_filename = report_name + '.html'
    porecov_params = {}
    tool_versions = {}
    tabledata = None
    col_formatters = {}
    col_descriptions = []


    def init(self, report_name):
        if report_name is not None:
            self.report_name = report_name
            self.output_filename = report_name + '.html'
        log(f'Created report object: {self.report_name}')


    def add_col_formatter(self, colname, colformatter):
        assert colname not in self.col_formatters, f'Duplicate column formatter: {colname}'
        self.col_formatters[colname] = colformatter


    def add_param(self, param_name, param_value):
        assert param_name not in self.porecov_params, f'Duplicate parameter: {param_name}'
        self.porecov_params[param_name] = param_value
        log(f'Added porecov param: {param_name}: {param_value}')


    def add_time_param(self):
        self.add_param('Report created', f'{time.strftime("%Y-%m-%d %H:%M:%S %Z", self.report_time)}')


    def add_version_param(self, porecov_version):
        warning_msg = 'Warning: Not an official release version of poreCov. Use parameter \'-r\' to specify a release version.'
        revision, commitID, scriptID = porecov_version.split(':')
        if revision != 'null':
            self.add_param('poreCov version', revision)
        else:
            if commitID != 'null':
                self.add_param('poreCov version', commitID + ' (git commitID) - ' + warning_msg)
            else:
                self.add_param('poreCov version', scriptID + ' (local scriptID) - ' + warning_msg)


    def parse_version_config(self, version_config_file):
        version_dict = {}
        try:
            with open(version_config_file) as infh:
                for line in infh:
                    lt = line.strip().replace(' ','')
                    if lt.startswith('withLabel:'):
                        tname, tinfo = lt.split(':',1)[1].split('{',1)
                        tcontainer = tinfo.split("'")[1] if "'" in tinfo else tinfo.split('"')[1]
                        tversion = tcontainer.split(':')[1].split('-')[0].lstrip('v')

                        assert tname not in version_dict, f'Duplicate tool in version config: {tname}'
                        version_dict[tname] = tversion
        except:
            print(version_dict)
            error(f'Failed to parse version config file: {version_config_file}')
        self.tool_versions = version_dict
        log('Parsed version config file.')


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


    def write_html_table(self, filehandle):
        filehandle.write(self.tabledata.to_html(classes=['tablestyle'], escape=False, formatters=self.col_formatters, float_format=lambda f: f'{f:.2f}'))


    def write_html_report(self):
        '''Write the html report to a file'''

        htmlheader = '''<!DOCTYPE html><html><head>
        <title>''' + self.report_name + '''</title>
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

        with open(self.output_filename, 'w') as outfh:
            outfh.write(htmlheader)
            outfh.write('<h1 class="header" id="main-header">poreCov Summary Report</h1>\n')

            # general params
            outfh.write('<h2 class="header" id="params-header">Run information</h2>\n')
            for param, value in self.porecov_params.items():
                outfh.write(param + ': \t' + value + '<br>\n')

            # results table
            outfh.write('<h2 class="header" id="table-header">Sample results</h2>\n')
            self.write_html_table(outfh)
            self.write_column_descriptions(outfh)

            outfh.write(htmlfooter)
        log(f'Wrote report to {self.output_filename}.')


    ### functions to add columns

    def add_pangolin_results(self, pangolin_results):
        log(f'Adding Pangolin results ...')
        # column names:
        # taxon,lineage,probability,pangoLEARN_version,status,note
        res_data = pd.read_csv(pangolin_results, index_col=0)
        self.check_and_init_tabledata(res_data.index)

        res_data['lineage_prob'] = [f'<b>{l}</b> ({p:.2f})' for l,p in zip(res_data['lineage'], res_data['probability'])]
        self.tabledata['Lineage (probability)'] = res_data['lineage_prob']

        self.add_col_description(f'Lineage and probability were determined with <a href="https://cov-lineages.org/pangolin.html">pangolin</a> (v{self.tool_versions["pangolin"]}).')
            

    def add_president_results(self, president_results):
        log(f'Adding President results ...')
        # column names:
        # ACGT Nucleotide identity	ACGT Nucleotide identity (ignoring Ns)	ACGT Nucleotide identity (ignoring non-ACGTNs)
        #	Date	Matches	Mismatches	N_bases	Ngap	acgt_bases	file_in_query	file_in_refiupac_bases	length_query
        #	length_reference	non_upac_bases	qc_all_valid	qc_post_align_pass_threshold	qc_post_aligned	qc_post_aligned_all_valid
        # 	qc_valid_length	qc_valid_nucleotides	qc_valid_number_n	query_description	query_index	query_name	reference_name
        res_data = pd.read_csv(president_results, index_col=24, sep='\t')
        self.check_and_init_tabledata(res_data.index)

        res_data['identity_mismatches'] = [f'{i*100:.2f} ({m})' for i, m in zip(res_data['ACGT Nucleotide identity'], res_data['Mismatches'])]
        self.tabledata['% Nucleotide identity (mismatches)'] = res_data['identity_mismatches']

        res_data['percentN_numN'] = [f'{nn/lr*100:.2f} ({nn})' for lr, nn in zip(res_data['length_reference'], res_data['N_bases'])]
        self.tabledata['%Ns (#Ns)'] = res_data['percentN_numN']

        self.add_col_description(f'Nucleotide identity, mismatches and Ns were determined with <a href="https://gitlab.com/RKIBioinformaticsPipelines/president">PRESIDENT</a> (v{self.tool_versions["president"]}).')

    
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

        def spike_markup(field):
            muts = field.split(', ')
            mumuts = []
            for mut in muts:
                if mut.split(':')[0] == 'S':
                    mumuts.append('<font color="#C54747"><b>' + mut + '</b></font>')
                else:
                    mumuts.append(mut)
            return ', '.join(mumuts)

        self.add_col_formatter('Mutations', spike_markup)
        self.add_col_description(f'Clade and mutations were determined with <a href="https://clades.nextstrain.org/">Nextclade</a> (v{self.tool_versions["nextclade"]}).')


    def add_kraken2_results(self, kraken2_results):
        log(f'Adding Kraken2 results ...')
        # column names:
        #sample,num_sarscov2,num_human
        res_data = pd.read_csv(kraken2_results, index_col=0)
        self.check_and_init_tabledata(res_data.index)

        res_data['total_reads'] = res_data['num_sarscov2'] + res_data['num_human']

        self.tabledata['%reads SARS-CoV-2'] = res_data['num_sarscov2'] / res_data['total_reads'] * 100.
        self.tabledata['%reads human'] = res_data['num_human'] / res_data['total_reads'] * 100.

        self.add_col_description(f'Read classification was determined with <a href="https://ccb.jhu.edu/software/kraken2/">Kraken2</a> (v{self.tool_versions["kraken2"]}).')


###

if __name__ == '__main__':

    log('Started summary_report.py ...')

    parser = argparse.ArgumentParser(description='Generate a summary report for multiple samples run with poreCov')
    parser.add_argument("-v", "--version_config", help="version config", required=True)
    parser.add_argument("--porecov_version", help="porecov_version", required=True)
    parser.add_argument("-p", "--pangolin_results", help="pangolin results")
    parser.add_argument("-n", "--nextclade_results", help="nextclade results")
    parser.add_argument("-q", "--president_results", help="president results")
    parser.add_argument("-k", "--kraken2_results", help="kraken2 results")
    args = parser.parse_args()


    # build report
    report = SummaryReport()

    # params
    report.parse_version_config(args.version_config)

    report.add_version_param(args.porecov_version)


    report.add_time_param()


    # results table, this determines the order of columns
    if args.kraken2_results:
        report.add_kraken2_results(args.kraken2_results)
    if args.pangolin_results:
        report.add_pangolin_results(args.pangolin_results)
    if args.president_results:
        report.add_president_results(args.president_results)
    if args.nextclade_results:
        report.add_nextclade_results(args.nextclade_results)

    
    report.write_html_report()


