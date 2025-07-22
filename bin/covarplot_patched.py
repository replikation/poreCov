#!/usr/bin/env python3
import os
import sys
import argparse
import gzip
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2021



    ----------------------------------------------------------------------------
    version 0.0.0 - initial
    version 0.0.1 - save plots and --log options added
    version 0.0.2 - variant depth calculations fixed


    TODO:
        - Colourblind options with --cb [1, 2, 3]
        - named amplicons from file or from bed
        - identify low coverage regions and highlight

    ----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2021 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def print_verbose(message):
    '''verbose printing'''
    sys.stderr.write('info: %s\n' % message)

def main():
    '''
    main func directing which plots to make
    '''
    VERSION = "0.0.3"
    NAME = "interArtic plots"
    parser = MyParser(
        description="Plots for interArtic",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # main options for base level checks and version output
    parser.add_argument("--version", action='version', version="{} version: {}".format(NAME, VERSION),
                        help="Prints version")

    parser.add_argument("-v", "--vcf_file",
                        help="full path to vcf file")


    parser.add_argument("-d1", "--depth_file_1",
                        help="full path to depth file 1")
    parser.add_argument("-d2", "--depth_file_2",
                        help="full path to depth file 2")

    parser.add_argument("-b", "--bed",
                        help="full path to scheme bed file")

    parser.add_argument("--show", action="store_true",
                        help="Show plot rather than saving it")
    parser.add_argument("-s", "--save",
                        help="Save path")
    parser.add_argument("-l", "--log", action="store_true",
                        help="y-axis log scale")


    args = parser.parse_args()
    print_verbose("arg list: {}".format(args))

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not args.show and not args.save:
        sys.stderr.write("No output method selected")
        parser.print_help(sys.stderr)
        sys.exit(1)

    bed_1, bed_2 = get_bed(args)

    if args.vcf_file and not args.depth_file_1:
        vcfx_snv, vcfy_snv, vcfx_id, vcfy_id = vcf_pipeline(args)
        plot(args, bed_1, bed_2, vcfx_snv=vcfx_snv, vcfy_snv=vcfy_snv, vcfx_id=vcfx_id, vcfy_id=vcfy_id)
    elif args.depth_file_1 and args.depth_file_2 and not args.vcf_file:
        both_covx, both_covy, cov1x, cov1y, cov2x, cov2y = cov_pipeline(args)
        plot(args, bed_1, bed_2, both_covx=both_covx, both_covy=both_covy,
        cov1x=cov1x, cov1y=cov1y, cov2x=cov2x, cov2y=cov2y)
    elif args.vcf_file and args.depth_file_1 and args.depth_file_2:
        vcfx_snv, vcfy_snv, vcfx_id, vcfy_id = vcf_pipeline(args)
        both_covx, both_covy, cov1x, cov1y, cov2x, cov2y = cov_pipeline(args)
        plot(args, bed_1, bed_2, vcfx_snv=vcfx_snv, vcfy_snv=vcfy_snv,
        vcfx_id=vcfx_id, vcfy_id=vcfy_id, both_covx=both_covx, both_covy=both_covy,
        cov1x=cov1x, cov1y=cov1y, cov2x=cov2x, cov2y=cov2y)
    else:
        sys.stderr.write("Command unknown: {}".format(args.command))
        parser.print_help(sys.stderr)
        sys.exit(1)

def get_bed(args):
    """
    create 2 lists of cordinates for overlapping amplicons
    """
    tmp_1 = []
    tmp_2 = []
    bed_1 = []
    bed_2 = []
    with open(args.bed, 'r') as f:
        for l in f:
            l = l.strip('\n')
            l = l.strip('\t')
            l = l.split('\t')
            #print(l)
            if "alt" in l[3]:
                continue
            if l[4][-1] == '1':
                tmp_1.append(int(l[1]))
                tmp_1.append(int(l[2]))
            elif l[4][-1] == '2':
                tmp_2.append(int(l[1]))
                tmp_2.append(int(l[2]))
            else:
                sys.stderr.write("bed format unknown: {}\n, please contact developers\n".format(l[-1]))

    tmp_1.sort()
    tmp_2.sort()

    for i in range(0,len(tmp_1)-3+1,4):
        bed_1.append((tmp_1[i], tmp_1[i+3]))
    for i in range(0,len(tmp_2)-3+1,4):
        bed_2.append((tmp_2[i], tmp_2[i+3]))


    return np.array(bed_1), np.array(bed_2)



def vcf_pipeline(args):
    """
    Plot vcf histogram
    x-axis = genome
    y-axis = depth
    each bar represents a variant, anotated
    Add genes/orfs/annoation to bottom (or maybe even top?) of plot for context
    """
    data = []
    x_snv = []
    y_snv = []
    x_id = []
    y_id = []
    header = []
    row = {}
    clair3 = False
    with gzip.open(args.vcf_file, 'rt') as f:
        for line_number, l in enumerate(f):
            
            # check if VCF is variant called with clair and we have to patch the INFO and FORMAT fields
            if line_number == 2 and 'Clair3' in l:
                clair3 = True
                print(f"Detected Clair3 as variant caller of VCF '{args.vcf_file}', patching variant depth compatibility.")

            if l[:2] == "##":
                continue
            if l[0] == "#":
                l = l[1:].strip('\n')
                l = l.split('\t')
                header = l
                continue

            l = l.strip('\n')
            l = l.split('\t')
            row = dict(zip(header, l))
            # how to calculate variant depths
            # medaka: depth = AC
            # nanopolish: depth = BaseCalledReadsWithVariant
            var_depth = -1

            if not clair3:
                for i in row["INFO"].split(";"):
                    # sys.stderr.write("vcf_INFO: {}\n".format(i))
                    if len(i) > 0:
                        name, result = i.split("=")
                        if name == "AC":
                            var_depth = int(result.split(",")[-1])
                        elif name == "BaseCalledReadsWithVariant":
                            var_depth = int(result)
            else:
                # for whatever reason clair3 puts the read depth into the FORMAT field and its value into SAMPLE:
                # ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">
                field_names  = row['FORMAT'].split(':')
                field_values = row['SAMPLE'].split(':')
                
                # check field is available and name is correct
                if len(field_names) >=2 and field_names[2] == 'DP' and len(field_values) >= 2:
                    var_depth = row['SAMPLE'].split(':')[2]

            if var_depth == -1:
                sys.stderr.write("var_depth not set for position: {}\n".format(row["POS"]))
            else:
                if len(row["REF"]) == len(row["ALT"]):
                    x_snv.append(int(row["POS"]))
                    y_snv.append(var_depth)
                else:
                    x_id.append(int(row["POS"]))
                    y_id.append(var_depth)

    if x_snv:
        nx_snv = np.array(x_snv)
        ny_snv = np.array(y_snv)
    if x_id:
        nx_id = np.array(x_id)
        ny_id = np.array(y_id)
    if x_snv and x_id:
        return nx_snv, ny_snv, nx_id, ny_id
    if x_snv and not x_id:
        return nx_snv, ny_snv, None, None
    if not x_snv and x_id:
        return None, None, nx_id, ny_id
    return None, None, None, None


def cov_pipeline(args):
    """
    Plot coverage histogram
    x-axis = genome
    y-axis = depth
    each bar is a bin of read depth
    """
    dp1x = []
    dp1y = []
    dp2x = []
    dp2y = []
    bothy = []

    with open(args.depth_file_1, 'r') as f:
        for l in f:
            l = l.strip("\n")
            l = l.split("\t")
            dp1x.append(int(l[2]))
            dp1y.append(int(l[3]))

    with open(args.depth_file_2, 'r') as f:
        for l in f:
            l = l.strip("\n")
            l = l.split("\t")
            dp2x.append(int(l[2]))
            dp2y.append(int(l[3]))

    for i in range(len(dp1y)):
        k = dp1y[i] + dp2y[i]
        bothy.append(k)

    bothy_trimmed = np.array([i for i in bothy])
    bothx_trimmed = np.array([i for i in dp1x])

    return np.array(bothx_trimmed), np.array(bothy_trimmed), np.array(dp1x), np.array(dp1y), np.array(dp2x), np.array(dp2y)

    # plt.fill_between(dp1x_trimmed, bothy_trimmed, color="skyblue", alpha=0.4)
    # plt.show()


def plot(args, bed_1, bed_2, vcfx_snv=None, vcfy_snv=None, vcfx_id=None, vcfy_id=None, both_covx=None, both_covy=None, cov1x=None, cov1y=None, cov2x=None, cov2y=None):
    """
    Plot everything separate or at once
    """
    if not args.show:
        plt.switch_backend('Agg')
    save_path = ""
    sample_name = ""
    fig, ax = plt.subplots()
    if args.vcf_file and not args.depth_file_1:
        # ax.stem(vcfx, vcfy, markerfmt= ' ')
        save_path = "/".join(args.vcf_file.split("/")[:-1])
        sample_name = args.vcf_file.split("/")[-1].split(".")[0]
        if vcfx_snv is not None:
            plt.vlines(vcfx_snv, 0, vcfy_snv, colors='darkorange')
        if vcfx_id is not None:
            plt.vlines(vcfx_id, 0, vcfy_id, colors='darkblue')
            plt.suptitle("Variants", y=0.95, fontsize=20)
    elif args.depth_file_1 and args.depth_file_2 and not args.vcf_file:
        save_path = "/".join(args.depth_file_1.split("/")[:-1])
        sample_name = args.depth_file_1.split("/")[-1].split(".")[0]
        if both_covx is not None:
            plt.fill_between(cov1x, cov1y, color="skyblue", alpha=0.6 , label="pool 1 coverage")
            plt.fill_between(cov2x, cov2y, color="violet", alpha=0.6, label="pool 2 coverage")
            plt.plot(both_covx, both_covy, color="darkgrey", alpha=0.6, linestyle="dotted", label="Combined coverage")
            plt.suptitle("Coverage", y=0.95, fontsize=20)
    elif args.vcf_file and args.depth_file_1 and args.depth_file_2:
        # ax.stem(vcfx, vcfy, markerfmt= ' ')
        save_path = "/".join(args.vcf_file.split("/")[:-1])
        sample_name = args.vcf_file.split("/")[-1].split(".")[0]
        if vcfx_snv is not None:
            plt.vlines(vcfx_snv, 0, vcfy_snv, colors='darkorange', label="SNV")
        if vcfx_id is not None:
            plt.vlines(vcfx_id, 0, vcfy_id, colors='darkblue', label="InDel")
        if both_covx is not None:
            plt.fill_between(cov1x, cov1y, color="skyblue", alpha=0.6, label="pool 1 coverage")
            plt.fill_between(cov2x, cov2y, color="violet", alpha=0.6, label="pool 2 coverage")
            plt.plot(both_covx, both_covy, color="darkgrey", alpha=0.6, linestyle="dotted", label="Combined coverage")
            # TODO: Make suptitle and title actually line up
            plt.suptitle("Variants and Coverage", y=0.95, fontsize=20)
    else:
        sys.stderr.write("this really shouldn't error. args: {}".format(args))


    #print(bed_1)
    # add amplicon to legend only once
    for i, j in bed_1[:1]:
        ax.add_patch(plt.Rectangle((i,-20),j-i, 15,facecolor='silver',
                              clip_on=False,linewidth = 1, label="amplicon"))
    for i, j in bed_1[1:]:
        ax.add_patch(plt.Rectangle((i,-20),j-i, 15,facecolor='silver',
                              clip_on=False,linewidth = 1))

    for i, j in bed_2:
        ax.add_patch(plt.Rectangle((i,-40),j-i, 15,facecolor='silver',
                              clip_on=False, linewidth = 1))

    plt.axhline(y=20, color='grey', linestyle='--', label="20x coverage min")

    if args.log:
        plt.yscale("log")
    plt.xlabel("Genome position", fontsize=20)
    plt.ylabel("Depth", fontsize=20)
    plt.tick_params(labelsize=10)
    plt.title(sample_name, fontsize=15)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    wi, hi = fig.get_size_inches()
    w = 1200
    h = 600
    fig.set_size_inches(hi*(w/h), hi)

    plt.tight_layout()
    if args.show:
        plt.show()
    elif args.save:
        save_file = args.save.rstrip("/") + "/" + sample_name + ".CoVarPlot.png"
        plt.savefig(save_file, dpi=1000/hi)
    else:
        save_file = save_path + "/" + sample_name + ".CoVarPlot.png"
        plt.savefig(save_file, dpi=1000/hi)



if __name__ == '__main__':
    main()
