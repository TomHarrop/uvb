#!/usr/bin/env python3

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

# run the Phytomine query
service = Service('https://phytozome.jgi.doe.gov/phytomine/service')
query = service.new_query('Gene')
query.add_view('name', 'primaryIdentifier', 'secondaryIdentifier',
               'organism.shortName', 'briefDescription')
query.add_constraint('Gene.name', 'ONE OF', gene_names)

# print the results into a tempfile
tmp = tempfile.mkstemp(suffix=".txt", text=True)[1]
sep = '\t'

with open(tmp, 'w') as outfile:
    outfile.write('name\tprimaryIdentifier\tsecondaryIdentifier\t'
                  'organism.shortName\tbriefDescription\n')
    for row in query.rows():
        line = (str(row['name']) + sep +
                str(row['primaryIdentifier']) + sep +
                str(row['secondaryIdentifier']) + sep +
                str(row['organism.shortName']) + sep +
                str(row['briefDescription']) + '\n')
        outfile.write(line)

# return filename for R to catch
print(tmp)

