#!/usr/bin/env python3
'''
Generate a summary report for multiple samples
./summary_report.py
'''
# SK

import sys
import time
import base64
import argparse
import pandas as pd


### functions

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\nAborting.\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')


### main

class SummaryReport():

    report_time = time.localtime()
    report_name = f'poreCov_summary_report_{time.strftime("%Y-%m-%d--%H-%M-%S", report_time)}'
    output_filename = report_name + '.html'
    porecov_params = {}
    tool_versions = {}
    scorpio_version = None
    scorpio_constellations_version = None
    variants_table = None
    pangolin_version = None
    pangolindata_version = None
    nextclade_version = None
    nextcladedata_version = None
    tabledata = None
    tabledataraw = None
    col_formatters = {}
    col_descriptions = []
    coverage_plots_b64 = []
    coverage_plots_filetype = []
    sample_QC_status = None
    sample_QC_info = {}
    control_string_patterns = ['control', 'negative']
    frameshift_warning = False
    samples_table = None


    # colors
    color_spike_markup = '#8a006d'
    color_good_green = '#046907'
    color_warn_orange = '#ac7800'
    color_error_red = '#a50500'


    def init(self, report_name):
        if report_name is not None:
            self.report_name = report_name
            self.output_filename = report_name + '.html'
        log(f'Created report object: {self.report_name}')


    def add_column(self, column_name, pandas_series):
        self.tabledata[column_name] = pandas_series

    def add_column_raw(self, column_name, pandas_series):
        self.tabledataraw[column_name] = pandas_series


    def add_col_formatter(self, colname, colformatter):
        assert colname not in self.col_formatters, f'Duplicate column formatter: {colname}'
        self.col_formatters[colname] = colformatter


    def add_param(self, param_name, param_value):
        assert param_name not in self.porecov_params, f'Duplicate parameter: {param_name}'
        self.porecov_params[param_name] = param_value
        log(f'Added porecov param: {param_name}: {param_value}')


    def add_QC_info(self, info_name, info_value):
        assert info_name not in self.sample_QC_info, f'Duplicate QC info: {info_name}'
        self.sample_QC_info[info_name] = info_value
        log(f'Added QC info: {info_name}: {info_value}')


    def add_time_param(self):
        self.add_param('Report created', f'{time.strftime("%Y-%m-%d %H:%M:%S %Z", self.report_time)}')


    def add_poreCov_version_param(self, porecov_version):
        pc_param = '<a href="https://github.com/replikation/poreCov"><b>poreCov</b></a> version'
        warning_msg = 'Warning: Not an official release version of poreCov. Use parameter \'-r\' to specify a release version.'
        revision, commitID, scriptID = porecov_version.split(':')
        if revision != 'null':
            self.add_param(pc_param, revision)
        else:
            if commitID != 'null':
                self.add_param(pc_param, commitID + ' (git commitID) - ' + warning_msg)
            else:
                self.add_param(pc_param, scriptID + ' (local scriptID) - ' + warning_msg)


    def add_pangolin_version_param(self):
        if self.pangolindata_version is None:
            error('add_pangolin_version_param() called before pangolin version was set')
        warning_msg = f' - <font color="{self.color_error_red}"><b>Warning</b>: A rather old version of pangolin-data was used ({self.pangolindata_version}). Use parameter \'--update\' to force the use of the most recent Pangolin container!</font>'
        
        pl_param = f'<a href="https://cov-lineages.org/resources/pangolin/requirements.html"><b>pangolin-data</b></a> version'
        pl_val = f'{self.pangolindata_version}'
        v1, v2, *rest = self.pangolindata_version.split('.')
        if int(v1) <= 1 and int(v2) <= 2:
            pl_val += warning_msg

        self.add_param(pl_param, pl_val)


    def add_nextclade_version_param(self):
        if self.nextcladedata_version is None:
            error('add_nextclade_version_param() called before parse_nextclade_version()')
        warning_msg = f' - <font color="{self.color_error_red}"><b>Warning</b>: A rather old version of Nextclade data was used ({self.nextcladedata_version}). Use parameter \'--update\' to force the use of the most recent Nextclade container!</font>'
        
        nc_param = f'<a href="https://clades.nextstrain.org/"><b>Nextclade</b></a> data version'
        nc_val = f'{self.nextcladedata_version}'

        year, month, day = self.nextcladedata_version.split('-')
        if int(year) <= 2021 and int(month) <= 10:
            nc_val += warning_msg

        self.add_param(nc_param, nc_val)


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


    def parse_nextclade_version(self, nextclade_docker):
        # e.g. nanozoo/nextclade:1.3.0--2021-06-25
        self.nextclade_version, self.nextcladedata_version = nextclade_docker.split(':',1)[-1].split('--')

    # UNUSED
    # def validate_index(self, t_index):
    #     '''Assert that an index matches the tabledata index.'''
    #     if self.tabledata is None:
    #         error('Cannot validate_index when tabledata is None')
    #     if not sorted(self.tabledata.index) == sorted(t_index):
    #         error(f'Index mismatch: tabledata: {sorted(self.tabledata.index)} \n query: {sorted(t_index)}')


    def check_and_init_table_with_samples(self, samples):

        print(samples)
        if samples == 'deactivated':
            log('No sample list input.')
        else:
            log('Using samples input.')
            self.samples_table = pd.read_csv(samples, index_col='_id')

            s_list = list(self.samples_table.index)
            log(f'Samples: {s_list}')

            s_table = pd.DataFrame(index=s_list)
            self.force_index_dtype_string(s_table)
            self.check_and_init_tabledata(s_table.index)

            # add description column if present
            if 'Description' in self.samples_table.columns:
                self.add_column_raw('Description', self.samples_table['Description'])
                self.add_column('Description', self.samples_table['Description'])


    def check_and_init_tabledata(self, t_index):
        '''If tabledata is None, initialize it now. Then check if all new index values are in the existing table.
        Thus samples input or adding the results with the most samples (kraken2) first is required.'''
        if self.tabledata is None:
            self.tabledata = pd.DataFrame(index=sorted(t_index))
            self.tabledata.columns.name = 'Sample'
            self.force_index_dtype_string(self.tabledata)
            self.add_col_description(f'Missing values (\'<font color="{self.color_error_red}">n/a</font>\') denote cases where the respective program could not determine a result.')
            self.tabledataraw = self.tabledata.copy()
        else:
            for item in t_index:
                assert item in self.tabledata.index, f'Index not found in existing table: {item}. Available: {self.tabledata.index}'


    def add_col_description(self, desc):
        self.col_descriptions.append(f'{desc}')


    def write_column_descriptions(self, filehandle):
        for desc in self.col_descriptions:
            filehandle.write(f'{desc}<br>\n')


    def write_html_table(self, filehandle):
        filehandle.write(self.tabledata.to_html(classes=['tablestyle'], escape=False, bold_rows=False, \
            na_rep=f'<font color="{self.color_error_red}">n/a</font>', formatters=self.col_formatters, float_format=lambda f: f'{f:.2f}'))


    def write_html_report(self):
        '''Write the html report to a file'''

        htmlheader = '''<!DOCTYPE html><html><head>
        <title>''' + self.report_name + '''</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">

        <style>
        * {
            font-family:"Open Sans",freesans,sans-serif
        }

        .content {
        max-width: 1700px;
        margin: auto;
        }

        table.tablestyle {
        background-color: #FFFFFF;
        text-align: center;
        border-collapse: collapse;
        }
        table.tablestyle td, table.tablestyle th {
        border: 2px solid #8B8B8B;
        padding: 5px 5px;
        }
        table.tablestyle tbody td {
        font-size: 16px;
        color: #000000;
        }
        table.tablestyle tr:nth-child(even) {
        background: #E6F5FF;
        }
        table.tablestyle thead {
        background: #E6F5FF;
        white-space: nowrap;
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
        summary {
        border: 1px solid #aaa;
        border-radius: 4px;
        white-space: nowrap;
        padding: 1px 4px;
        background-color: #eee;
        }
        code {
        background-color: #eee;
        border-radius: 3px;
        font-family: courier, monospace;
        padding: 0 3px;
        }

        </style>

        <script>
        const expandElements = (shouldExpand, whatToExpand) => {
            let detailsElements = document.querySelectorAll("details[tag=" + whatToExpand + "]");
            
            detailsElements = [...detailsElements];

            if (shouldExpand) {
                detailsElements.map(item => item.setAttribute("open", shouldExpand));
            } else {
                detailsElements.map(item => item.removeAttribute("open"));
            }
        };
        </script>

        </head>

        <body>
        <div class="content">
        '''

        htmlfooter = '''
        </div>
        </body></html>
        '''

        with open(self.output_filename, 'w') as outfh:
            outfh.write(htmlheader)
            outfh.write('<h1 class="header" id="main-header">poreCov Summary Report</h1>\n')

            # general
            outfh.write('<h2 class="header" id="params-header">Run information</h2>\n')
            for info, value in self.sample_QC_info.items():
                outfh.write(value + '<br>\n')
            outfh.write('<br>\n')
            for param, value in self.porecov_params.items():
                outfh.write(param + ': ' + value + '<br>\n')

            # results table
            outfh.write('<h2 class="header" id="table-header">Sample results</h2>\n')
            self.write_html_table(outfh)
            self.write_column_descriptions(outfh)

            if self.coverage_plots_b64 != []:
                self.write_html_coverage_plot(outfh)

            outfh.write(htmlfooter)
        log(f'Wrote report to {self.output_filename}.')


    ### functions to add columns

    def force_index_dtype_string(self, dataframe):
        dataframe.index = dataframe.index.astype('string')
        

    def add_pangolin_results(self, pangolin_results):
        log(f'Adding Pangolin results ...')
        # column names used:
        # taxon,lineage,conflict,scorpio_call,scorpio_conflict
        res_data = pd.read_csv(pangolin_results, index_col='taxon', dtype={'taxon': str})
        self.force_index_dtype_string(res_data)
        self.check_and_init_tabledata(res_data.index)

        # pangolin and scorpio versions
        # column names used:
        # version,pangolin_version,scorpio_version,constellation_version
        assert res_data.shape[0] > 0
        self.pangolin_version = res_data.iloc[0]['pangolin_version']
        self.pangolindata_version = res_data.iloc[0]['version'].split('-',1)[-1].split('v',1)[-1]
        self.scorpio_version = res_data.iloc[0]['scorpio_version']
        self.scorpio_constellations_version = res_data.iloc[0]['constellation_version']


        # get data
        self.add_column_raw('pangolin_lineage', res_data['lineage'])
        self.add_column_raw('pangolin_conflict', res_data['conflict'])

        res_data['lineage_conflict'] = [f'<b>{l}</b><br>({p if pd.notnull(p) else "-"})' for l,p in zip(res_data['lineage'], res_data['conflict'])]
        colname = 'Lineage<br>(conflict)'
        
        # get VOC/VOI status if present
        if self.variants_table is not None:
            colname += f'<br>Variant type'
            for row in res_data.index:
                var_status = self.get_lineage_status(res_data.at[row, 'lineage'])
                res_data.at[row, 'variant_status'] = var_status

                # annotate pangolin column
                if var_status != '':
                    color = self.color_warn_orange
                    if var_status in ['VOI', 'VOC']:
                        color = self.color_error_red
                    res_data.at[row, 'lineage_conflict'] += f'<br><font color="{color}"><b>{var_status}</b></font>'

            self.add_column_raw('variant_status', res_data['variant_status'])
            self.add_col_description(f'Variant type (VOC, VOI, etc.) was determined from the <a href="{args.variants_table}">variants table</a> ' + \
                'of <a href="https://github.com/3dgiordano/SARS-CoV-2-Variants">SARS-CoV-2-Variants</a>.')


        self.add_column(colname, res_data['lineage_conflict'])
        if self.pangolin_version is None or self.pangolindata_version is None:
            error('No pangolin/pangolin-data versions were added before adding pangolin results.')
        self.add_col_description(f'Lineage and the corresponding tree resolution conflict measure were determined with ' + \
            f'<a href="https://cov-lineages.org/pangolin.html">Pangolin</a> (v{self.pangolin_version} using ' + \
            f'<a href="https://cov-lineages.org/resources/pangolin/requirements.html">pangolin-data</a> v{self.pangolindata_version}).')
        
        # Add scorpio info if any is present
        if res_data['scorpio_call'].notna().any():
            log(f'Found "scorpio_call" value(s), adding Scorpio results ...')
        
            self.add_column_raw('scorpio_constellation', res_data['scorpio_call'])
            self.add_column_raw('scorpio_conflict', res_data['scorpio_conflict'])

            res_data['scorpio_conflict'] = [f'<b>{l if pd.notna(l) else "-"}</b><br>({p if pd.notna(p) else "-"})' for l,p in zip(res_data['scorpio_call'], res_data['scorpio_conflict'])]

            self.add_column('Constellation<br>(conflict)', res_data['scorpio_conflict'])
            if self.scorpio_version is None or self.scorpio_constellations_version is None:
                error('No Scorpio/constellations versions were added before adding Pangolin results.')
            self.add_col_description(f'Constellation and the corresponding tree resolution conflict measure were determined with ' + \
                f'<a href="https://github.com/cov-lineages/scorpio">Scorpio</a> (v{self.scorpio_version} using ' + \
                f'<a href="https://cov-lineages.org/constellations.html">Constellations</a> version {self.scorpio_constellations_version}).')



    def add_president_results(self, president_results):
        log(f'Adding President results ...')

        res_data = pd.read_csv(president_results, index_col='query_name', sep='\t', dtype={'query_name': str})
        self.force_index_dtype_string(res_data)
        self.check_and_init_tabledata(res_data.index)


        # identity and mismatches
        self.add_column_raw('president_identity', res_data['ACGT Nucleotide identity'])
        self.add_column_raw('president_mismatches', res_data['Mismatches'])

        def identity_markup(ident, mismatches):
            color = self.color_good_green
            if ident < 99.:
                color = self.color_warn_orange
            if ident < 90.:
                # RKI rule
                color = self.color_error_red
            return  f'<font color="{color}">{ident:.2f}</font><br>(<font color="{color}">{int(mismatches)}</font>)'

        res_data['identity_mismatches'] = [identity_markup(i*100, m) if not pd.isnull(m) else m for i, m in zip(res_data['ACGT Nucleotide identity'], res_data['Mismatches'])]
        self.add_column('%ident.<br>(mism.)', res_data['identity_mismatches'])


        # percent and number Ns
        res_data['percN'] = res_data['N_bases']/res_data['length_query']*100
        self.add_column_raw('president_percentN', res_data['percN'])
        self.add_column_raw('president_numN', res_data['N_bases'])

        def percN_markup(nn, ql):
            color = self.color_good_green
            percN = nn/ql*100
            if percN > 1.:
                color = self.color_warn_orange
            if percN > 5.:
                # 5% RKI rule
                color = self.color_error_red
            return f'<font color="{color}">{percN:.2f}</font><br>(<font color="{color}">{nn}</font>)'

        res_data['percN_numN'] = [percN_markup(nn, ql) for nn, ql in zip(res_data['N_bases'], res_data['length_query'])]

        self.add_column('%Ns<br>(#Ns)', res_data['percN_numN'])
        self.add_col_description(f'Nucleotide identity, mismatches, Ns and QC pass status were determined with <a href="https://gitlab.com/RKIBioinformaticsPipelines/president">PRESIDENT</a> (v{self.tool_versions["president"]}).')


        # QC pass column & save QC status
        if self.sample_QC_status is not None:
            error('sample_QC_status is already set when running add_president_results')
        self.sample_QC_status = {}

        for sample in self.tabledata.index:
            if sample in res_data.index and res_data.loc[sample, 'qc_all_valid']:
                res_data.loc[sample, 'QC_pass'] = f'<font color="{self.color_good_green}"><b>YES</b></font>'
                res_data.loc[sample, 'QC_pass_raw'] = 'YES'
                self.sample_QC_status[sample] = 'pass'
            else:
                res_data.loc[sample, 'QC_pass'] = f'<font color="{self.color_error_red}"><b>NO</b></font>'
                res_data.loc[sample, 'QC_pass_raw'] = 'NO'
                self.sample_QC_status[sample] = 'fail'

        self.add_column_raw('president_QC_pass', res_data['QC_pass_raw'])
        self.add_column('QC<br>pass', res_data['QC_pass'])
        self.add_col_description(f'QC pass criteria are the RKI genome submission requirements: >= 90% identity to NC_045512.2, <= 5% Ns, etc. (<a href="https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Qualitaetskriterien.pdf?__blob=publicationFile">PDF</a> in German).')

    
    def add_nextclade_results(self, nextclade_results):
        log(f'Adding Nextclade results ...')

        res_data = pd.read_csv(nextclade_results, index_col='seqName', sep='\t', dtype={'seqName': str})
        self.force_index_dtype_string(res_data)
        self.check_and_init_tabledata(res_data.index)


        self.add_column_raw('nextclade_clade', res_data['clade'])
        self.add_column_raw('nextclade_mutations_nt', res_data['substitutions'])
        self.add_column_raw('nextclade_mutations', res_data['aaSubstitutions'])
        self.add_column_raw('nextclade_deletions_nt', res_data['deletions'])
        self.add_column_raw('nextclade_deletions', res_data['aaDeletions'])
        self.add_column_raw('nextclade_insertions_nt', res_data['insertions'])
        self.add_column_raw('nextclade_insertions', res_data['aaInsertions'])
        self.add_column_raw('nextclade_frameshifts', res_data['frameShifts'])

        # private mutation information
        self.add_column_raw('nextclade_privateNucMutations_reversion', res_data["privateNucMutations.reversionSubstitutions"])
        self.add_column_raw('nextclade_privateNucMutations_labeled', res_data["privateNucMutations.labeledSubstitutions"])
        self.add_column_raw('nextclade_privateNucMutations_unlabeled', res_data["privateNucMutations.unlabeledSubstitutions"])
        self.add_column_raw('nextclade_privateNucMutations_qc_status', res_data["qc.privateMutations.status"])

        # N information
        self.add_column_raw('nextclade_missing', res_data["missing"])

        res_data['mutations_formatted'] = [m.replace(',', ', ') if type(m) == str else '-' for m in res_data['aaSubstitutions']]
        res_data['deletions_formatted'] = [m.replace(',', ', ') if type(m) == str else '-' for m in res_data['aaDeletions']]
        res_data['insertions_formatted'] = [m.replace(',', ', ') if type(m) == str else '-' for m in res_data['aaInsertions']]
        res_data['frameshifts_formatted'] = [m.replace(',', ', ') if type(m) == str else '-' for m in res_data['frameShifts']]

        if (res_data['frameshifts_formatted'] != '-').any():
            self.frameshift_warning = True


        def get_button_html(tag):
            return f'<button onClick="expandElements(true, \'{tag}\')">Expand</button><button onClick="expandElements(false, \'{tag}\')">Collapse</button>'

        self.add_column('Clade', res_data['clade'])
        muts_colname = f'Mutations<br>(<font color="{self.color_spike_markup}"><b>on spike</b></font>)<br>' + get_button_html('mutations')
        dels_colname = f'Deletions<br>(<font color="{self.color_spike_markup}"><b>on spike</b></font>)<br>' + get_button_html('deletions')
        inss_colname = f'Insertions<br>(<font color="{self.color_spike_markup}"><b>on spike</b></font>)<br>' + get_button_html('insertions')
        frms_colname = f'Frameshifts<br>(<font color="{self.color_spike_markup}"><b>on spike</b></font>)<br>' + get_button_html('frameshifts')
        self.add_column(muts_colname, res_data['mutations_formatted'])
        self.add_column(dels_colname, res_data['deletions_formatted'])
        self.add_column(inss_colname, res_data['insertions_formatted'])
        self.add_column(frms_colname, res_data['frameshifts_formatted'])

        def clade_markup(field):
            return f'<b>{field}</b>'

        def spike_markup(field):
            muts = field.split(', ')
            mumuts = []
            for mut in muts:
                first = mut.split(':')[0]
                if first == 'S' or (first.isnumeric() and 21563 <= int(first) <= 25385):
                    mumuts.append(f'<font color="{self.color_spike_markup}"><b>' + mut + '</b></font>')
                else:
                    mumuts.append(mut)
            return ', '.join(mumuts)

        def spike_markup_with_toggle(field, tag):
            if field == '-':
                return field
            else:
                n_mutations = int(len(field.split(',')))
                return f'<details tag={tag}><summary>Number found: <b>{n_mutations}</b></summary>\n{spike_markup(field)}</details>'

        def get_markup_with_toggle_and_tag(tag):
            return lambda x: spike_markup_with_toggle(x, tag)


        self.add_col_formatter('Clade', clade_markup)
        self.add_col_formatter(muts_colname, get_markup_with_toggle_and_tag('mutations'))
        self.add_col_formatter(dels_colname, get_markup_with_toggle_and_tag('deletions'))
        self.add_col_formatter(inss_colname, get_markup_with_toggle_and_tag('insertions'))
        self.add_col_formatter(frms_colname, get_markup_with_toggle_and_tag('frameshifts'))

        if self.nextclade_version is None or self.nextcladedata_version is None:
            error('No nextclade/nextcladedata versions were added before adding nextclade results.')
        self.add_col_description(f'Clade, mutations, deletions, insertions and frameshifts were determined with ' + \
            f'<a href="https://clades.nextstrain.org/">Nextclade</a> (v{self.nextclade_version} using nextclade data release {self.nextcladedata_version}).')
    
        

    def add_kraken2_results(self, kraken2_results):
        log(f'Adding Kraken2 results ...')
        # column names:
        # sample,num_unclassified,num_sarscov2,num_human
        res_data = pd.read_csv(kraken2_results, index_col='sample', dtype={'sample': str})
        self.force_index_dtype_string(res_data)
        self.check_and_init_tabledata(res_data.index)

        def readable_si_units(number):
            e3steps = 0
            units = ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
            while number > 1e3 and e3steps < 8:
                e3steps += 1
                number /= 1e3
            return f'{number:.1f}{units[e3steps]}' if e3steps > 0 else f'{number:.0f}'

        # colors
        def sars_markup(value):
            color = self.color_good_green
            if value < 95.:
                color = self.color_warn_orange
            if value < 90.:
                color = self.color_error_red
            return  f'<font color="{color}">{value:.2f}</font>'

        def human_markup(value):
            color = self.color_good_green
            if value > 1.:
                color = self.color_warn_orange
            if value > 5.:
                color = self.color_error_red
            return  f'<font color="{color}">{value:.2f}</font>'

        def unclass_markup(value):
            color = self.color_good_green
            if value > 5.:
                color = self.color_warn_orange
            if value > 10.:
                color = self.color_error_red
            return  f'<font color="{color}">{value:.2f}</font>'

        res_data['total_reads'] = res_data['num_unclassified'] + res_data['num_sarscov2'] + res_data['num_human']
        perc_sarscov_colname = '%reads<br>SARS-CoV-2<br>(#reads)'
        perc_human_colname = '%reads<br>human<br>(#reads)'
        perc_unclass_colname = '%reads<br>unclass.<br>(#reads)'

        self.add_column_raw('kraken2_numreads_sarscov2', res_data['num_sarscov2'])
        self.add_column_raw('kraken2_numreads_human', res_data['num_human'])
        self.add_column_raw('kraken2_numreads_unclassified', res_data['num_unclassified'])
        self.add_column_raw('kraken2_numreads_total', res_data['total_reads'])

        res_data['n_sars'] = [f"{sars_markup(n_sars/n_total*100.)}<br>({readable_si_units(n_sars)})" \
            for n_sars, n_total in zip(res_data['num_sarscov2'], res_data['total_reads'])]
        res_data['n_human'] = [f"{human_markup(n_human/n_total*100.)}<br>({readable_si_units(n_human)})" \
            for n_human, n_total in zip(res_data['num_human'], res_data['total_reads'])]
        res_data['n_unclass'] = [f"{unclass_markup(n_unclass/n_total*100.)}<br>({readable_si_units(n_unclass)})" \
            for n_unclass, n_total in zip(res_data['num_unclassified'], res_data['total_reads'])]

        self.add_column(perc_sarscov_colname, res_data['n_sars'])
        self.add_column(perc_human_colname, res_data['n_human'])
        self.add_column(perc_unclass_colname, res_data['n_unclass'])

        self.add_col_description(f'Read classification was determined against a database containing only SARS-CoV-2 and human with <a href="https://ccb.jhu.edu/software/kraken2/">Kraken2</a> (v{self.tool_versions["kraken2"]}).')


    def add_coverage_plots(self, coverage_plots):
        for coverage_plot in sorted(coverage_plots.split(',')):
            data = open(coverage_plot, 'rb').read() # read bytes from file
            data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
            self.coverage_plots_b64.append(data_base64.decode())
            ftype = coverage_plot.rsplit('.',1)[-1]
            if ftype not in ('png', 'jpg', 'jpeg', 'gif'):
                log(f'WARNING: {coverage_plot} does not look like a good image filetype.')
            self.coverage_plots_filetype.append(ftype)


    def write_html_coverage_plot(self, filehandle):
        if self.coverage_plots_b64 is []:
            error('No coverage plot was added beforehand.')
        filehandle.write(f'''<h2>Coverage plots</h2>
        Coverage of all samples against the SARS-CoV-2 reference genome (NC_045512.2) determined with <a href="https://github.com/lh3/minimap2">minimap2</a> (v{self.tool_versions["minimap2"]}) and <a href="https://github.com/RaverJay/fastcov">fastcov</a> (v{self.tool_versions["fastcov"]}).<br>
        ''')
        for plot, ftype in zip(self.coverage_plots_b64, self.coverage_plots_filetype):
            filehandle.write(
            f'''<img src="data:image/{ftype};base64,{plot}"><br>
            ''')


    def check_if_control(self, sample_name):
        for pattern in self.control_string_patterns:
            if pattern in sample_name:
                return True
        return False


    def add_QC_status_info(self):
        if self.sample_QC_status is None:
            error('sample_QC_status was not set before calling add_QC_status_info().')

        n_realsamples = 0
        n_controls = 0
        n_passrealsamples = 0
        n_passcontrols = 0
        for sample, status in self.sample_QC_status.items():
            if self.check_if_control(sample):
                n_controls += 1
                if status == 'pass':
                    n_passcontrols += 1
            else:
                n_realsamples += 1
                if status == 'pass':
                    n_passrealsamples += 1

        # add status info
        if n_passrealsamples > 0:
            self.add_QC_info('Passed samples', f'<font color="{self.color_good_green}"><b>{n_passrealsamples} / {n_realsamples} of samples passed QC criteria.</b></font>')
        if n_passrealsamples < n_realsamples:
            self.add_QC_info('Failed samples', f'<font color="{self.color_error_red}"><b>{n_realsamples-n_passrealsamples} / {n_realsamples} of samples failed QC criteria.</b></font>')
        if n_controls > 0:
            if n_passcontrols < n_controls:
                self.add_QC_info('Negative controls', f'<font color="{self.color_good_green}"><b>{n_controls-n_passcontrols} / {n_controls} of control samples correctly did not produce an assembly that passed QC criteria.</b></font>')
            if n_passcontrols > 0:
                self.add_QC_info('Bad controls', f'<font color="{self.color_error_red}"><b>{n_passcontrols} / {n_controls} of control samples wrongly produced an assembly that passed QC criteria.</b></font>')
        else:
            self.add_QC_info('Negative control', f'<font color="{self.color_warn_orange}"><b>No negative control samples were found by automatic detection.</b></font>')
        patterns = "'" + "', '".join(self.control_string_patterns) + "'"
        self.add_QC_info('Note', f'Note: samples are considered negative controls if their name contains certain keywords ({patterns}) - please check if these assignments were correct.')


        # warning if variant table is missing
        if self.variants_table is None:
            self.add_QC_info('Var warning', f'<font color="{self.color_warn_orange}">Warning: variants table was not provided - potential information on VOCs, VOIs, etc. is missing from the lineage column.</font>')

        # warning if variant table fallback was used
        if self.variants_table is not None and 'fallback' in args.variants_table:
            self.add_QC_info('fallback warning', f'<font color="{self.color_warn_orange}">Warning: variants table could not be downloaded and a fallback was used (<a href="{args.variants_table}">{args.variants_table}</a>). Please check if it is outdated! To fix, your setup need to be able to download from github.</font>')

        # frameshift warning message
        if self.frameshift_warning:
            self.add_QC_info('Frameshift warning', f'<font color="{self.color_error_red}"><b>WARNING:</b> There were frameshifts in one or more samples. ' + \
                'This is most likely an error due to incorrect consensus calling as frameshifts are biologically very unlikely. ' + \
                '<b>CAUTION:</b> This error can cause masking of downstream mutations, deletions and insertions in affected protein regions! ' + \
                'Basecalling with superior accuracy (\'sup\') models and/or using --nanopolish for consensus calling can fix these errors.</font>')


        # mark control samples
        def mark_controls(sample_name):
            if self.check_if_control(sample_name):
                return sample_name + f'<br>(<font color="{self.color_spike_markup}">considered control</font>)'
            else:
                return sample_name

        self.tabledata.index = [mark_controls(sn) for sn in self.tabledata.index]


    def write_table_output(self):
        if self.tabledataraw is None:
            error('Failed to write table data - no raw table data was added beforehand.')
        self.tabledataraw.to_excel(self.report_name + '_datatable.xlsx' , sheet_name='poreCov', index_label='sample')
        self.tabledataraw.to_csv(self.report_name + '_datatable.tsv', index_label='sample', sep='\t')


    def parse_variants_table(self, table_file):
        try:
            self.variants_table = pd.read_csv(table_file, header=0)
            log(f'Read in variants table: {table_file}')
        except:
            log(f'Could not read variants table: {table_file}')


    def get_lineage_status(self, lineage):
        if self.variants_table is None:
            error('get_lineage_status() used before parse_variants_table()')
        
        # label 	pango 	interest 	type
        # Alpha 	B.1.1.7 	WHO 	VOC
        # Alpha - Q.1 	Q.1 	WHO 	VOC
        # ...

        candidates = self.variants_table[self.variants_table['pango'] == lineage]['type']
        if len(candidates) > 1:
            log(f'Multiple hits for lineage {lineage}')
            return 'ambigious'
        if len(candidates) == 1:
            return candidates.iloc[0] if candidates.notnull().all() else ''

        # check if there is a .* mask

        # B.1.177 -> B.1.* ?
        mask = lineage.rsplit('.', 1)[0] + '.*'
        candidates = self.variants_table[self.variants_table['pango'] == mask]['type']
        if len(candidates) > 1:
            log(f'Multiple hits for mask {mask}')
            return 'ambigious'
        if len(candidates) == 1:
            return candidates.iloc[0] if candidates.notnull().all() else ''

        # B.1.177 -> B.1.177.* ?
        mask = lineage + '.*'
        candidates = self.variants_table[self.variants_table['pango'] == mask]['type']
        if len(candidates) > 1:
            log(f'Multiple hits for mask {mask}')
            return 'ambigious'
        if len(candidates) == 1:
            return candidates.iloc[0] if candidates.notnull().all() else ''

        return 'n/a'




### main

if __name__ == '__main__':

    log('Started summary_report.py ...')

    parser = argparse.ArgumentParser(description='Generate a summary report for multiple samples run with poreCov')
    parser.add_argument("-v", "--version_config", help="version config", required=True)  
    parser.add_argument("--variants_table", help="variants table with VOCs, VOIs etc.", required=True)    
    parser.add_argument("--porecov_version", help="porecov version", required=True)
    parser.add_argument("--guppy_used", help="guppy used")
    parser.add_argument("--guppy_model", help="guppy model")
    parser.add_argument("--medaka_model", help="medaka model")
    parser.add_argument("--nf_commandline", help="full nextflow command call", required=True)
    parser.add_argument("--nextclade_docker", help="nextclade/nextcladedata version", required=True)
    parser.add_argument("--primer", help="primer version")
    parser.add_argument("-p", "--pangolin_results", help="pangolin results")
    parser.add_argument("-n", "--nextclade_results", help="nextclade results")
    parser.add_argument("-q", "--president_results", help="president results")
    parser.add_argument("-k", "--kraken2_results", help="kraken2 results")
    parser.add_argument("-c", "--coverage_plots", help="coverage plots (comma separated)")
    parser.add_argument("-s", "--samples", help="sample ids (comma separated)")
    args = parser.parse_args()


    ### build report
    report = SummaryReport()
    report.parse_version_config(args.version_config)
    report.parse_variants_table(args.variants_table)
    report.parse_nextclade_version(args.nextclade_docker)


    # check for samples input
    if args.samples:
        report.check_and_init_table_with_samples(args.samples)


    # results table, this determines the order of columns
    if args.kraken2_results:
        report.add_kraken2_results(args.kraken2_results)
    if args.president_results:
        report.add_president_results(args.president_results)
    if args.pangolin_results:
        report.add_pangolin_results(args.pangolin_results)
    if args.nextclade_results:
        report.add_nextclade_results(args.nextclade_results)
    if args.coverage_plots:
        report.add_coverage_plots(args.coverage_plots)

    # metadata

    # total QC status
    report.add_QC_status_info()

    # params
    report.add_poreCov_version_param(args.porecov_version)
    report.add_pangolin_version_param()
    report.add_nextclade_version_param()

    # check run type
    if args.primer:
        
        # infer run type from guppy usage
        report.add_param('<br>Run type', "Genome reconstruction and classification from raw sequencing data " + \
            ("(fast5)" if args.guppy_used == 'true' else "(fastq)") + (' utilizing Nanopolish (fast5) information' if '--nanopolish' in args.nf_commandline else ''))
        report.add_param('<a href="https://artic.network/ncov-2019">ARTIC</a> version', report.tool_versions['artic'])
        report.add_param('ARTIC primer version', args.primer)
        # add guppy/medaka model if used
        if args.guppy_model and args.guppy_used == 'true':
            report.add_param('Guppy model', args.guppy_model)
        if args.medaka_model:
            report.add_param('Medaka model', args.medaka_model)
    else:
        report.add_param('<br>Run type', "Genome classification from sequences (fasta)")
    report.add_param('poreCov command', '<code>' + args.nf_commandline + '</code>')
    report.add_time_param()

    report.write_html_report()
    report.write_table_output()

