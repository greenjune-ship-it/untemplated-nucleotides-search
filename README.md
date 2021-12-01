# test-pipeline

## Requirements
* bowtie 1.2.3
* samtools 1.10
* graphviz 2.43.0
* cutadapt 2.8 with Python 3.8.10
* seqtk 1.3
* BBMap and default-jre
* GPL Ghostscript 9.50
* weblogo 3.7.8

## How to run the pipeline on example dataset

```
cd example
snakemake -j 4 --snakefile ../Snakefile
```

For dry run without execution you can use
```
snakemake -j 4 --dry-run
```

## General recommendation about running the pipeline on you own dataset
I recommend you to keep the following directories stucture
```
|-analysis_name
  |
  |-data
  |  |
  |  |-reference*
  |  |-reads*
  |
  |-configs
  | |
  | |-config.yaml*
  |
  |-results
    |
    |--iter_0
    |  |
    |  |-alignments
    |  |-fastq
    |  |-original
    |
    |--iter_1
    |--iter_n
```

FoFolders ```data/reference``` and ```data/reads``` should be preexisting and specified in config file.


### Required files
```configs/config.yaml``` must containg all necessary information about pipeline parameters (see example/configs/congig.yaml) 

```data/reference``` must contain reference genome in fasta format

```data/reads``` must contain reads in fastq format

Other folders will appear during the execution of pipeline.

### Outputs
```results/iter_*/alignments``` contains all mapped results, including mapped and unmapped alignments in bam format

```results/iter_*/fastq``` contains fastq files extracted from bam

```results/iter_*/original``` contains original non-trimmed reads in fasta fromat and their sequence logos for each length specified in config file

## Create pipeline graph

```
snakemake --dag | dot -Tsvg > dag.svg
```

<p align="center">
  <img src="./example/dag.svg">
</p>

```
snakemake --rulegraph | dot -Tsvg > rulegraph.svg
```

<p align="center">
  <img src="./example/rulegraph.svg">
</p>
