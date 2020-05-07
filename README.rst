===========================================================
snakePipes
===========================================================

.. image:: https://readthedocs.org/projects/snakepipes/badge/?version=latest
    :target: http://snakepipes.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.org/maxplanck-ie/snakepipes.svg?branch=develop
    :target: https://travis-ci.org/maxplanck-ie/snakepipes
    :alt: Build Staus

.. image:: https://zenodo.org/badge/54579435.svg
    :target: https://zenodo.org/badge/latestdoi/54579435
    :alt: Citation


snakePipes are flexible and powerful workflows built using `snakemake <snakemake.readthedocs.io>`__ that simplify the analysis of NGS data.

.. image:: ./docs/content/images/snakePipes.png
   :scale: 20 %
   :height: 100px
   :width: 100 px
   :align: right

Workflows available
--------------------

- DNA-mapping*
- ChIP-seq*
- RNA-seq*
- ATAC-seq*
- scRNA-seq
- Hi-C
- Whole Genome Bisulfite Seq/WGBS

**(*Also available in "allele-specific" mode)**

We have made modifications to the DNA-mapping and ChIP-seq workflows, in order to make them compatible with the current practices in the Zwart lab. These include mapping single-end reads with BWA, marking duplicates with picard, filtering out reads with a MAPQ score below 20, and estimating fragment size with phantompeakquals for MACS2 peakcalling of single-end reads. The others workflows remain intact from the original pipeline, and their functionality should not be affected.

Installation
-------------

Ensure conda is properly installed by running:

``conda --version``

If the command does not work, you may need to initialize conda in your environment with:

``conda init``

You may need to log out and log back in for changes to take effect.

Once you have conda working, clone this repository into your desired location with:

``git clone https://github.com/csijcs/snakepipes.git``

Change directory into the snakepipes folder with:

``cd snakepipes``

Then run the following:

``conda env create --file snakepipes.yaml``

This will create a new conda environment called "snakepipes" into which snakePipes is installed. You will now need to create the conda environments needed by the various workflows.

First activate the snakepipes environment with:

``conda activate snakepipes``

Then run the build script with:

``sh build.sh``

You now need to create the various environments required for the pipeline by running:

``snakePipes createEnvs --condaDir ~/.conda/envs``

Before running any anlyses, you will need to create indices. These only need to be created once, but each genome build (i.e. hg19, hg38, mm10, etc.) will need their own indices.  These can be constructed with the following command:

``createIndices --genomeURL <path/URL to your genome fasta> --gtfURL <path/url to genes.gtf> --local -o <output_dir> <name>``. 

The necessary files/links can be obtained from https://www.gencodegenes.org/

For example, to create the required indicies for hg19 the command would be:

``createIndices --genomeURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz --gtfURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.transcripts.fa.gz --local -o /PATH/TO/OUTPUT/DIRECTORY/hg19 hg19``

to create the required indicies for hg38 the command would be:

``createIndices --genomeURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz --gtfURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz --local -o /PATH/TO/OUTPUT/DIRECTORY/hg38 hg38``

to create the required indicies for mm10 the command would be:

``createIndices --genomeURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz --gtfURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz --local -o /PATH/TO/OUTPUT/DIRECTORY/mm10 mm10``

You will need to supply your own /PATH/TO/OUTPUT/DIRECTORY/ above (i.e. the location where you want the genome indices stored). 

Once indices are created, you are ready to proceed to the pipelines.

DNA-mapping

For the DNA-mapping pipeline, the minimum required command is:

``DNA-mapping -i /INPUT/DIR -o /OUTPUT/DIR --local genome_build`` 

The default mapping program is Bowtie2. To use BWA, copy the above bwa_mapping.yaml to the directory you are running the pipeline from. For example, for mapping to with BWA to hg19, first put all .fastq.gz files into a folder named FASTQ. Then run the following command:

``DNA-mapping -i /PATH/TO/FASTQ -o /PATH/TO/OUTPUT/DIRECTORY --configfile bwa_mapping.yaml --local -j 10 --mapq 20 --trim --trim_prg cutadapt --fastqc hg19``

Here, -i specifies the input folder contaning the .fastq.gz files, -o is the output directory, --local runs on the local server and not on a cluster, -j specifies the number of threads, --trim tells the pipeline to trim the reads, --trim_prg tells the pipeline the program used to trim the reads, --fastqc tell it to run fastqc analysis, and finally hg19 specifies the genome.

ChIP-seq

The ChIP-seq pipline is designed to take the ouput directly from the DNA-mapping pipeline. The only additional file you will need is a sample_config.yaml file, telling the progrom your sample names, the control for each sample, and whether they to look for broad peaks (i.e. histone marks) or narrow peaks (i.e. transcription factors). See the example sample_config.yaml file above.

If you have run the DNA-mapping pipeline first, then simply run:

``ChIP-seq -d /PATH/TO/DNA-mapping/OUTPUT --local -j 10 --single-end hg19 sample_config.yaml``

Here -d should be the directory with the output of the DNA-mapping pipeline, and it will also direct the output of the ChIP-seq pipeline there. If your samples are not single end then remove the --single-end flag. Also modify the genome_build (i.e. hg19) to suit your purposes).

If you have not run the DNA-mapping pipeline first, then you can still run the pipeline directly from BAM files. In this case, put all of your .bam files into a folder called "bams" (or whatever you want). You will also need the from_bam.yaml file from above in the working directory. Additional parameters (such as fragment length) can also be modified in this file. Then run:

``ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bams --configfile from_bam.yaml --local -j 10 --single-end hg19 sample_config.yaml``

There will be various folder outputs, including some QC, but the peak files will be in the MACS2 folder. In the future we will likely implement some additional QC measures, such as cross-correlation ("phantom peaks"), and possibly add modules for DiffBind and other downstream analysis. For now this will get the reads mapped and peaks called effectively.

The other modules have remained untouched and should work according to the original pipeline.

Documentation
--------------

For detailed documentation on setup and usage, please visit our `read the docs page <https://snakepipes.readthedocs.io/en/latest/>`__.


Citation
-------------

If you adopt/run snakePipes for your analysis, cite it as follows :

Bhardwaj V, Heyne S, Sikora K, Rabbani L, Rauer M, Kilpert F, et al. **snakePipes enable flexible, scalable and integrative epigenomic analysis.** bioRxiv. 2018. p. 407312. `doi:10.1101/407312 <https://www.biorxiv.org/content/early/2018/09/04/407312>`__


Note
-------------

SnakePipes are under active development. We appreciate your help in improving it further. Please use issues to the GitHub repository for feature requests or bug reports.
