#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
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

# read the taxon ids
taxon_ids = {}
with open('data/phytozome_species_ids.txt') as csv_file:
    t_ids = csv.reader(csv_file, delimiter=',')
    for row in t_ids:
        taxon_ids[row[0]] = row[1]

# read the gene names
with open(gene_name_list, 'r') as f:
    gene_names = [line.strip() for line in f]

# run the Phytomine query
views = ['primaryIdentifier', 'organism.shortName',
         'goAnnotation.ontologyTerm.identifier',
         'goAnnotation.ontologyTerm.name',
         'goAnnotation.ontologyTerm.namespace']
service = Service('https://phytozome.jgi.doe.gov/phytomine/service')
query = service.new_query('Gene')
query.add_views(views)
query.add_constraint('Gene.name', 'ONE OF', gene_names, code='A')
query.add_constraint('organism.taxonId', 'ONE OF',
                     list(taxon_ids.values()), code='B')
query.set_logic("A and B")

# write results
tmp = tempfile.mkstemp(suffix=".txt", text=True)[1]
sep = "\t"
with open(tmp, 'w') as outfile:
    # write header
    outfile.write(sep.join(views) + '\n')
    for row in query.rows():
        # write a line from row with the result for each view
        outfile.write(sep.join(list(str(row[x]) for x in views)) + '\n')

# return filename for R to catch
print(tmp)
