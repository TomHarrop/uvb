#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile
import argparse
from intermine.webservice import Service

# read the gene IDs from --input
parser = argparse.ArgumentParser(
    description='Fetch gene information from Phytomine')
parser.add_argument('--input', '-i', help='List of Gene.names (one per line)',
                    type=str, required=True, dest='gene_name_list')
args = parser.parse_args()
gene_name_list = args.gene_name_list

with open(gene_name_list, 'r') as f:
    gene_names = [line.strip() for line in f]

# set up phytomine query
service = Service('https://phytozome.jgi.doe.gov/phytomine/service')
views = ['primaryIdentifier', 'homolog.gene2.primaryIdentifier',
         'homolog.gene2.organism.shortName', 'homolog.relationship']
query = service.new_query('Gene')
query.add_views(views)
query.add_constraint('Gene.name', 'ONE OF', gene_names)
query.add_constraint('homolog.gene2.organism.id', 'ONE OF',
                     ['171000000', '256000000', '276000000', '328000000',
                      '204000000', '241000000'], code='B')
query.set_logic('A and B')

# print the results into a tempfile as tsv
tmp = tempfile.mkstemp(suffix=".txt", text=True)[1]
sep = '\t'
with open(tmp, 'w') as outfile:
    # write header
    outfile.write(sep.join(views) + '\n')
    for row in query.rows():
        # write a line from row with the result for each view
        outfile.write(sep.join(list(str(row[x]) for x in views)) + '\n')

# return filename for R to catch
print(tmp)
