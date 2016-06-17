
---
link-citations: true
...

# Transcriptional response of land plants to environmental doses of UV-B radiation

Thomas Harrop^\*,1^, Laura Monforte López^†^, Swarit Jasial^\*,2^, Roman Ulm^‡^, Javier Martinez Abaigar^†^ and Ales Pecinka^\*,§^

^\*^ Department of Plant Breeding and Genetics, Max Planck Institute for Plant Breeding Research, 50829 Cologne, Germany

^†^ Department of xxxxxxxxxxxxx

^‡^ Department of xxxxxxx

^§^ For correspondance. Ph.: +49 221 5062 465; E-mail: pecinka@mpipz.mpg.de

^1^ Present address: UMR DIADE, Institut de Recherche pour le Développement, 34394 Montpellier, France

^2^ Present address: 

Sequence data from this article are available from the National Centre for Biotechnology Information Sequence Read Archive under accession number **SRP0000000**.

Running title: 

Keywords: plants, evolution, UV-B, transcriptomics, XXXXXXX

## Abstract

TF

## Introduction

* Role of UV-B in plant development
* Exposure of land plants to UV-B
* Choice of plant species for study

## Results and discussion

We used whole-genome expression analysis to investigate the effect of environmental doses of UV-B radiation on transcription in different plant species. Plants were exposed to 15 min of broad spectrum (BS) or high wavelength (HW) UV-B light and allowed to recover for two hours in photosynthetically active radiation (PAR) conditions. RNA was isolated from plant tissues and used to generate Illumina TruSeq polyA+ libraries for RNA-sequencing. Two separate control treatments were also performed, where plants either remained under PAR or were covered with 320 nm high-pass filters during UV-B exposure. Two biological replicates were performed for each treatment and control, and a single sequencing library was produced for each replicate. At least 20 M 50 base, single-ended reads were generated for each library. After adaptor trimming and quality filtering, an average of 85% of reads from each library uniquely mapped to the appropriate reference genome, and 38 of the 56 libraries produced greater than 20 M uniquely-mapping reads within annotated genes (**supplementary data**). Only reads from *O. sativa* and *S. moellendorffii* libraries mapped at rates lower than 75%, with a corresponding decrease in the number of reads mapped to genes (Figure S1). This may be explained by gaps in the reference sequence or annotation, or technical factors including inefficent removal of rRNA from sequencing libraries and contamination with sequences from other species.

### Transcriptional response to UV-B radiation

We tested for differential expression between control libraries and either BS or HW UV-B at a Log~2~-fold change (L~2~FC) level of 0.58 with a false discovery rate of 0.1. DESeq2 detected six differentially expressed genes (DEGs) in *A. thaliana*, eight in each of *C. rheinhardtii* and *S. lycopersicum* and one in *O. sativa* and *S. moellendorffii* [**supplementary data**; @Love:2014do]. Although most of the DEGs are not functionally characterized, there are some exeptions. The *A. thaliana* *FLAVANONE 3-HYDROXYLASE* gene, which **does what** in regulation and synthesis of flavonoids (**no TAIR access... look for references**), was induced by BS and HW UV-B, as was a putative *O. sativa* phenylalanine ammonia-lyase, *LOC_Os04g43800*. Two other chlorophyll A-B binding protein genes, which have been previously shown (**check refs**) to be inducible by light, *EARLY LIGHT-INDUCIBLE PROTEIN 2* (*A. thaliana*) and *Cre17.g740950* (*C. rheinhardtii*), were also induced.

* Possible reasons for the low number of genes detected using this method:
    - combined controls (what happens if I only compare against PAR). Why is it important to use both controls?
    - genuinely low response to this level of UV-B (compare to number of DEGs in microarray paper)
    - low number of replicates reduces power to detect DEGs

Because the power to detect differentially expressed genes is limited when using a low number of biological replicates (**ref?**), a clustering approach was also used to recover further genes of interest. For each species, the normalised expression values were transformed using the variance-stabilising transformation (VST) included with DESeq2 [@Love:2014do; @Anders:2010fu], filtered to remove lowly-expressed genes (maximum VST-expression value < 8), and sorted to retrieve the top 10% of genes by variance. The expression values for these genes were standardised and subjected to a soft clustering step using Mfuzz 2.26.0 [@Kumar:2007uw]. To visualise the clusters, the expression of each gene was plotted in two dimensions by subtracting the geometric mean of its VST-expression values in the four control libraries from the geometric mean of its VST-expression values in the libraries generated from plants exposed either to broad spectrum or high-wavelength UV-B, and the point for each gene was coloured according to the cluster it was assigned to (figure 1).

* Clustering responses:
    - aquatic species, Sp and Cr, don't respond a lot
    - globally, most genes that responded to both treatments did so in the same direction. Only At and Os had genes that were induced by one treatment and downregulated by the other

The results of this analysis reveal several patterns in the response of different species to UV-B. The two aquatic species, *C. reinhardtii* and *S. polyrhiza*, had the smallest response to the two treatments in terms of magnitude, suggesting that under the conditions used in these experiments they are not as sensitive to UV-B as the other species at the level of gene regulation. *S. moelendorffi* had the largest response, and clusters corresponding to genes that were up-regulated by both treatments or specifically up- or down-regulated by broad-spectrum UV-B are evident. In two species, tomato and *P. patens*, many genes were induced by both treatments, suggesting that there is a large overlap in response to the two conditions used in this experiment. However, there is also a cluster of tomato genes that were down-regulated by high-wavelength UV-B but up-regulated by broad-spectrum UV-B, suggesting that the transcriptional machinery was able to differentiate between the treatments. This response was not observed in P. patens, although there are genes that are specifically down-regulated by broad-spectrum UV-B. Most other species had a more modular response to UV-B, with evidence of clusters of genes that responded differently to the two treatments.

To investigate potential biological meanings of the expression clusters, the list of genes from each cluster was analysed using the DAVID data-mining environment [@Huang:2009gk]. This analysis was only possible for species that are present in the DAVID database, namely *A. thaliana*, *O. sativa* and *P. patens*. The terms found in each cluster were compared within each of these species using hierarchical clustering of distance matrices based on the enriched terms for each cluster (see methods). 

In *A. thaliana*, genes that were more strongly up-regulated by broad-spectrum UV-B than by high-wavelength UV-B (clusters 1 and 6) were enriched for genes annotated with the term ‘phenypropanoid biosynthetic process’ (figure 2). Genes up-regulated by high-wavelength UV-B (clusters 2, 3 and 4) were enriched for ‘hormone-mediated signalling’, ‘response to jasmonic acid stimulus’, ‘response to gibberellin stimulus’ and ‘regulation of transcription’. This suggests that hormone production and signalling is an important part of the response to UV-B and underlines the significance of the differential expression of genes involved in flavonoid synthesis.

The results in *O. sativa* were less clear, although the cluster of genes down-regulated after exposure to broad-spectrum UV-B but moderately up-regulated by high-wavelength UV-B contain genes with annotations relating to photosynthesis, and the genes annotated with the term ‘dicarboxylic acid metabolic process’ also contain a moderately enriched cluster of genes annotated with the term ‘biosynthesis of phenylpropanoids’ (figure 3). Similarly, one cluster of genes upregulated by both treatments in P. patens (cluster 6) contains the term ‘monosaccharide metabolic process’, but these genes are also annotated with terms relating to hormone synthesis incuding ‘biosynthesis of phenylpropanoids’ (figure 4). DAVID analysis is fundamentally dependent on classification of protein functions, which may not be perfect [@Schnoes:2009gb], particularly in less well-annotated species such as rice. However, these results suggest that hormone signalling, especially by phenypropanoids such as flavonoids, are an important component of the plant’s response to UV-B in several species.

### Inter-specific comparison

* Establishment and description of homology groups
* Venn diagramming at inter-specific level (or some matrix or heat map?)
* Focus on flavonol pathway?
* More detailed analysis of selected pathway?

## Materials and Methods

### Plant material

The plant material is listed in Table X. For Arabidopsis, rice and tomato, seeds were first surface sterilized and then grown on sterile ½ Murashige-Skoog media. S. polyrhiza was cultivated in liquid XXXX media (recipe from Petra). S. moellendorffii plants were grown on a peat soil mixed with stone in closed containers to ensure high air humidity. P. patens was maintained on xxxxxxxx media as a protonema culture and gametophores were induced by addition of XXXXXXX. All species were cultivated in a growth chamber at 21°C under long day (16 h light : 8 h dark) conditions. PAR instensity???. Transcriptomic analysis were performed on the first and the second pair of true leaves of Arabidopsis, rice and tomato, entire plants of S. polyrhiza, vegetative branches of S. moellendorffii and gametophores of P. patens. (add vegetative phylloids of M. polymorpha?).

### UV treatments

Plants were covered with specific glass cut off filters WG 295 nm, 305 nm and 320 nm (Schott) and UV irradiated using four Philips TL40 W / 12RS SLV UV fluorescent tubes (fluence rate at 312 nm: 4.35 W / cm²) for 15 min. The UV treatment was started 6 h after the onset of the chamber day. Control plants were treated in a UV free chamber section. After the treatment, plants were recovered for 2 h at UV-free conditions.

### Nucleic acids isolation, cDNA synthesis and RNA-sequencing

RNA was extracted using the RNeasy kit (Qiagen), with additional on column DNase I (Roche) treatment. cDNA for qPCR experiments was synthesized from 1 μg RNA per sample with Revert Aid H-Minus First Strand  cDNA synthesis kit using the oligo-d(T) primer (Thermo Scientific). The purity of cDNA was monitored by PCR with an intron-spanning primer pair.

RNA sequencing was performed with two biological replicates per experimental point. The libraries were prepared from 1 μg total RNA with RNA integrity number >7.8 (Bioanalyzer; Agilent) using TruSeq RNA kit (Illumina) and sequenced as 100 bp single-end reads on HiSeq2500 (Illumina). Reads were trimmed and low quality reads filtered with FAST-X tools (http://hannonlab.cshl.edu/fastx_toolkit/) using custom made scripts. This yielded an average of 15 million high quality reads per library. The reads were mapped to the corresponding reference genome assembly and annotation (Table X) using the Tophat2 (REFERENCE) with default settings. The coverage of individual genes was retrieved with the Qualimap from the set of uniquely mapped reads and significance (adjusted P-value = <0.05) of transcriptional changes estimated with the DEseq package (REFERENCE) in R.

### Quantitative PCR

The RT-qPCR was performed using 1 μl cDNA per 10 µl reaction with SensiMix (PeqLab) kit on an CFX384 (Biorad) instrument. Fold-changes were calculated relative to mock treated controls using the standard curve method. Primers used in this study are provided in Supplementary Table X.

### Homology searches

Here, I definitely need your help!

### Data availability

Sequence data from this article are available from the National Centre for Biotechnology Information Sequence Read Archive under accession number **SRP0000000**.

Other MM parts?

## Acknowledgments

We thank following colleagues for donations of plant material: B. Reiss for P. patens, K. Yokawa for O. sativa, S. Rensing for S. moellendorffii and K. Appenroth for S. polyrhiza. We are grateful for technical assistance of B. Eilts, P. Pecinkova and R. Gentges and critical reading of the manuscript by M. Koornneef. This work was supported by the funding from Max Planck Society to A.P. and T.H. and the XXXX grant XXXXX to L.M.L.

## Conflicts of interest

The authors declare no conflicts of interest.

## Tables

Table 1. Significantly differentially expressed genes (adjusted p-value < 0.1 & |fold change| > 1.5).

## Figure legends

Figure 1. Change in expression after exposure to high-wavelength or broad-spectrum UV-B. For each gene, the geometric mean of the transformed expression values in the control libraries was subtracted from the geometric mean of transformed expression values after exposure to either broad-spectrum (x-axis) or high-wavelength (y-axis) UV-B. Genes are coloured according to their membership to clusters that were defined by soft clustering of expression values. Note that the cluster number is not intended to correspond between species. The responses to the treatments vary between different species.

Figure 2. Results from DAVID analysis of gene expression clusters in A. thaliana. Following the use of functional annotation clustering with DAVID, the annotation groups were compared between expression clusters to pick representative terms from each annotation group (y-axis; see methods). The tiles of the heatmap are coloured by the enrichment score of the annotation group, and labeled with the number of genes within each annotation cluster. Note that even blue tiles represent enrichment, but the enrichment score in these cases is less than 1.3, which corresponds to a p-value of 0.05 {Huang:2009gk}.

Figure 3. Results from DAVID analysis of gene expression clusters in O. sativa. See methods and figure 2.

Figure 4. Results from DAVID analysis of gene expression clusters in P. patens. See methods and figure 2.

## Literature cited


