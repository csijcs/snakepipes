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

I have currently only tested our modifications for the DNA-mapping and ChIP workflows. The others remain intact from the original pipeline, and their functionality should not be affected.

Installation
-------------

You will need python > 3.5, so first check with:

``python --version``

Snakepipes uses conda for installation and dependency resolution, so you will need to `install conda <https://conda.io/docs/user-guide/install/index.html>`__ first.

Ensure conda is properly installed by running:

``conda --version``

Afterward, simply run the following:

``conda create -n snakePipes -c mpi-ie -c bioconda -c conda-forge snakePipes``

This will create a new conda environment called "snakePipes" into which snakePipes is installed. You will then need to create the conda environments needed by the various workflows.

First run:
``conda activate snakePipes`` to activate the appropriate conda environment.

Then run:
``snakePipes createEnvs`` to create the various environments and register GATK.

Indices and annotations needed to run the workflows can be created by a simple command :

``createIndices --genomeURL <path/URL to your genome fasta> --gtfURL <path/url to genes.gtf> --local -o <output_dir/genome_build>``. 

For example, to create the required indicies for hg19 the command would be:

``createIndices --genomeURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz --gtfURL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz --local -o /PATH/TO/OUTPUT/DIRECTORY/hg19``

Indices only need to be created once. You are now ready to proceed to the piplines.

DNA-mapping
For the DNA-mapping pipeline, the minimum required command is:

``DNA-mapping -i /INPUT/DIR -o /OUTPUT/DIR --local genome_build`` 
(The --local flag is not expressly required, however if you are running on a local server and not a cluster then it will be necessary).

For example, for mapping to hg19, first put all .fastq.gz files into a folder named FASTQ. Then run the following command:

``DNA-mapping -i /PATH/TO/FASTQ -o /PATH/TO/OUTPUT/DIRECTORY --local -j 10 --mapq 20 --trim --trim_prg cutadapt --fastqc hg19``

Here, -i specifies the input folder contaning the .fastq.gz files, -o is the output directory, --local runs on the local server and not on a cluster, -j specifies the number of threads, --trim tells the pipeline to trim the reads, --trim_prg tells the pipeline the program used to trim the reads, --fastqc tell it to run fastqc analysis, and finally hg19 specifies the genome.

ChIP-seq
The ChIP-seq pipline is designed to take the ouput directly from the DNA-mapping pipeline. The only additional file you will need is a sample_config.yaml file, telling the progrom your sample names, the control for each sample, and whether they to look for broad peaks (i.e. histone marks) or narrow peaks (i.e. transcription factors). See the example sample_config.yaml file above.

If you have run the DNA-mapping pipeline first, then simply run:

``ChIP-seq -d /PATH/TO/DNA-mapping/OUTPUT --local -j 10 --single-end hg19 sample_config.yaml``

Here -d should be the directory with the output of the DNA-mapping pipeline, and it will also direct the output of the ChIP-seq pipeline there. If your samples are not single end then remove the --single-end flag. Also modify the genome_build (i.e. hg19) to suit your purposes).

If you have not run the DNA-mapping pipeline first, then you can still run the pipeline directly from BAM files. In this case, put all of your .bam files into a folder called "bams" (or whatever you want) and run:

``ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bams --local -j 10 --single-end hg19 sample_config.yaml``

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
