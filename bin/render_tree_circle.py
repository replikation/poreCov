#!/usr/bin/env python3

import argparse
import toytree      
import toyplot
import numpy as np
import toyplot.pdf
import toyplot.svg
from io import StringIO

parser = argparse.ArgumentParser(
    description='''
    help for help
    ''')

parser.add_argument(
    '--tree', required=True,
    help='Tree in newick format')

parser.add_argument(
    '--highlight', required=True,
    help='which leaf to highlight')

parser.add_argument(
    '--format', required=True,
    help='which treeformat')

args = parser.parse_args()

newick = args.tree
highlighttip = args.highlight
formattree = args.format

# Creating tree

tree = toytree.tree(newick, tree_format = int(formattree))

# Tree Style

mystyle = {
    "layout": 'circular',
    "edge_type": 'p',
    "edge_style": {
        "stroke": "black",
        "stroke-width": 1,
    },
    "tip_labels_align": True,
    #"tip_labels_colors": "black",
    #"node_labels": False,
    "edge_align_style": {
        "stroke": "grey",
        "stroke-width": 0.5,
        "stroke-dasharray": "2"
    }, 
}

colorlist = ["#d6557c" if highlighttip in tip else "#5384a3" for tip in tree.get_tip_labels()]


canvas, axes = tree.draw(
    tip_labels_colors=colorlist,
    height=(tree.ntips * 8),
    width=(tree.ntips * 8),
    **mystyle,
    tip_labels_style={"font-size": "10px"},
    node_labels=None,
    node_sizes=5,
    node_colors='grey'
);

# stretch domain to fit long tip names
# axes.x.domain.max = 5

toyplot.svg.render(canvas, "tree_circle.svg")
toyplot.pdf.render(canvas, "tree_circle.pdf")