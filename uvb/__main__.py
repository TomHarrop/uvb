#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
#################
# UV-B pipeline #
#################

import csv
import functions
import ruffus
import os

######################
# PIPELINE FUNCTIONS #
######################


# download genomes from phytozome
def download_genome(output_files, species, genome_url, annotation_url,
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


# generate STAR indices
def generate_index(input_files, output_files, species):
    job_script = 'src/sh/generate_index'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = species + "_index"
    job_id = functions.submit_job(job_script, ntasks, cpus_per_task, job_name,
                                  extras=['-s', species])
    functions.print_job_submission(job_name, job_id)


# define reads
def define_reads(output_files, species):
    path_to_reads = "data/reads/" + species
    assert os.path.isdir(path_to_reads), ("Error: reads folder " +
                                          path_to_reads + " missing")
    read_files = os.listdir(path_to_reads)
    print("[ Using reads in folder " + path_to_reads + " ]")
    for fileName in read_files:
        qual_name = path_to_reads + '/' + fileName
        assert os.path.isfile(qual_name), ("Error: read file " + qual_name +
                                           " missing")
        print(qual_name)
    functions.touch(output_files)


# run cutadapt
def cutadapt(input_reads, output_reads, species):
    job_script = 'src/sh/cutadapt'
    ntasks = '4'
    cpus_per_task = '1'
    job_name = species + '_cutadapt'
    job_id = functions.submit_job(job_script, ntasks, cpus_per_task, job_name,
                                  extras=['-s', species])
    functions.print_job_submission(job_name, job_id)


# run star mapping
def star(input_files, output_files, species):
    job_script = 'src/sh/star'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = species + '_star'
    job_id = functions.submit_job(job_script, ntasks, cpus_per_task, job_name,
                                  extras=['-s', species])
    functions.print_job_submission(job_name, job_id)


# run deseq2
# def deseq(input_files, output_files):
#    pass

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
        main_pipeline.originate(
            name=species + "_genome",
            task_func=download_genome,
            output="data/genome/" + species + "/METADATA.csv",
            extras=[species, fasta_urls[species], annotation_urls[species],
                    jgi_logon, jgi_password])

    # generate a star genome for each species
    star_indices = main_pipeline.transform(
        task_func=generate_index,
        input=ruffus.output_from(list(species + "_genome" for species in
                                      fasta_urls.keys())),
        filter=ruffus.regex(r"data/genome/(.*)/METADATA.csv"),
        output=r"output/\1/star-index/METADATA.csv",
        extras=[r"\1"])

    # define the reads
    for species in fasta_urls.keys():
        main_pipeline.originate(name=species + "_reads",
                                task_func=define_reads,
                                output="ruffus/" + species + "_reads",
                                extras=[species])

    # first mapping step
    for species in fasta_urls.keys():
        main_pipeline.collate(
            name=species + "_mapped_reads",
            task_func=star,
            input=[[ruffus.output_from(species + "_reads"),
                    "output/" + species + "/star-index/METADATA.csv"]],
            filter=ruffus.formatter(),
            output=["output/{subdir[1][1]}/star/METADATA.csv"],
            extras=["{subdir[1][1]}"])\
            .follows(star_indices)

    # downstream somethingorother
#    deseq_results = main_pipeline.transform(
#                task_func=deseq,
#                input=ruffus.output_from(list(species + "_mapped_reads" for species in fasta_urls.keys())),
#                filter=ruffus.regex(r"output/(\w+)/star/METADATA.csv"),
#                output=[r"output/\1/deseq/METADATA.csv"],
#                extras=[r"\1"])

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

    # print the flowchart
    ruffus.pipeline_printout_graph("ruffus/flowchart.pdf", "pdf",
                                   pipeline_name="UV-B analysis pipeline")

if __name__ == "__main__":
    main()
