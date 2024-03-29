## Draft

https://docs.google.com/document/d/1VQN_7j0omnzKNc566jfe39stZVhXsfG2feNEkaJ4cFc/edit

## Goals

* Optimize GFF format definition and usability
* Detect methodology accuracy due to tools and some experimental step in the protocols.

## Data

### Tewari data

https://www.ncbi.nlm.nih.gov/pubmed/30010675

http://www.biorxiv.org/content/biorxiv/early/2017/05/17/113050.full.pdf

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=whipakmajrwprcv&acc=GSE94586

### Narry Kim data

https://academic.oup.com/nar/article/47/5/2630/5271499

### Carrie Wrigth data 

https://www.biorxiv.org/content/10.1101/445437v3

### DSRG data

Still to be published, another study to compare protocols using the mirxplor sample.

### Fratta data

Evaluation of methodologies for microRNA biomarker detection by next generation sequencing
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6161688/
issue: wrong ID number to download data. Contacted the author to get the data.


### To be added if we get the data

Systematic assessment of commercially available low-input miRNA library preparation kits
https://www.biorxiv.org/content/10.1101/702456v1.full
No data yet.

## Processed data

Trimming was done with the [smadann](https://github.com/miRTop/smadann) nexflow pipeline.

The following command was used for each study and type of data:

`nextflow run mirtop/smadann --csv totrim.csv -c ../../om-profile.config --outdir trimmed -qs 10`

Analysis was done with [bowtie]  + [mirtop] in a [snakemake] file located in each study and data type.

`snakemake -p -s run.snakefile`

### standard analysis of synthetic

Mirxplor reference was parsed to use only synthetic with an edit distance of 4 or more, and the alignments were filtered to keep only reads that mapped to those unique synthetic with a maximum of 4 changes. Code used for this is at [analysis](analysis/analysis.R) folder.

Data is available for anyone at [aws mirtop space](https://mirtop-tewari-data.s3.amazonaws.com/synthetic_2019_srr_mirgff1.2.rds). 

Currently contains: tewari, wrigth, kim and dsrg data.

### biological samples

For human data we use miRBase22 to map all sequences. Same filtering step were used here.

Data is available for anyone at [aws mirtop space](https://mirtop-tewari-data.s3.amazonaws.com/real_2019_srr_mirgff1.2.rds). 

## Tools

* bcbio smallRNA-seq pipeline + isomiRs - On charge Lorena Pantano
* isomiR-SEA - On charge Gianvito Urgese
* ChimiRa, miRge - On charge Marck Halushka
* sRNAbench - On charge Michael Hackenberg
* Prost - Thomas Desvignes
* miRGe - Marc K. Halushka
* (Add your tool here and person will do it)

## Questions to address

* Reproducibility of replicates
* Reproducibility of protocols
* Reproducibility of tools

## Results

Updated report can be found [here](results/results.md)

## Milestones:

### Set up

* [X] Select random public data
* [X] Run with all the tools listed above
* [X] Put data in common space
* [X] Adapt output tools to GFF format

### Random sample

Sample [SRR5756178](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5756178) is a whole blood small RNA-seq run from this manuscript https://academic.oup.com/nar/article/4080663 and is part of project PRJNA391912.  It has ~ 2.8 million reads, of which ~2.6 million are miRNAs.

### Synthetic data

Benchmark was done with synthetic isomiRs for one human miRNA, see [results](https://github.com/miRTop/incubator/tree/master/synthetic).

