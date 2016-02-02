#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
#################
# UV-B pipeline #
#################

import csv
import functions
import ruffus


######################
# PIPELINE FUNCTIONS #
######################

# download genomes from phytozome
def download_genome(outputFiles, species, genome_url, annotation_url,
                    jgi_logon, jgi_password):
    job_script = 'src/sh/download_genomes'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = species + "_dl"
    jgi_prefix = 'http://genome.jgi.doe.gov'
    job_id = functions.submit_job(job_script, ntasks, cpus_per_task, job_name,
                                  extras=['-e', jgi_logon, '-p', jgi_password,
                                          '-s', species,
                                          '-g', jgi_prefix + genome_url,
                                          '-a', jgi_prefix + annotation_url])
    functions.print_job_submission(job_name, job_id)


#########################
# PIPELINE CONSTRUCTION #
#########################

def main():

    # prepare the ruffus pipeline
    main_pipeline = ruffus.Pipeline.pipelines["main"]

    # catch jgi logon and password from cli
    parser = ruffus.cmdline.get_argparse(description='UV-B analysis pipeline.')
    parser.add_argument('--email', '-e',
                        help='Logon email address for JGI',
                        type=str,
                        dest='jgi_logon')
    parser.add_argument('--password', '-p',
                        help='JGI password',
                        type=str,
                        dest='jgi_password')
    options = parser.parse_args()
    jgi_logon = options.jgi_logon
    jgi_password = options.jgi_password

    # need a dictionary of species to genome URL and species to gff.
    # supply this in a text file
    fasta_urls = {}
    annotation_urls = {}
    with open('data/genomeUrls.txt') as tsv:
        genome_urls = csv.reader(tsv, delimiter='\t')
        next(genome_urls, None)
        for row in genome_urls:
            fasta_urls[row[0]] = row[1]
            annotation_urls[row[0]] = row[2]

    # iterate over fasta_urls keys to call download script
    for species in fasta_urls.keys():
        genomes = main_pipeline.originate(
            name=species + "_genome", task_func=download_genome,
            output="data/genome/" + species + "/METADATA.csv",
            extras=[species, fasta_urls[species], annotation_urls[species],
                    jgi_logon, jgi_password])

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

    # print the flowchart
    ruffus.pipeline_printout_graph("ruffus/flowchart.pdf", "pdf",
                                   pipeline_name="UV-B analysis pipeline")

if __name__ == "__main__":
    main()
