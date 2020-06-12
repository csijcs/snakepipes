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


This repository is a modification of the snakePipes workflows forked from `https://github.com/maxplanck-ie/snakepipes`

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

Begin by logging into the server and entering the terminal environment.

First initialize conda in your environment with:

``conda init``

Ensure conda is properly installed by running:

``conda --version``

Find the path to the conda your are using with:

``which conda``

For harris server (most if not all users), your conda path should be:

``/opt/anaconda/bin/conda``

If you are on darwin, your conda path should be:

``/opt/anaconda3/bin/conda``

Once you have conda working and the proper path, clone this repository into your desired location with:

``git clone https://github.com/csijcs/snakepipes.git``

Change directory into the snakepipes folder with:

``cd snakepipes``

Run the following using your pull conda path (harris shown for example):

``/opt/anaconda/bin/conda env create --file snakepipes.yaml``

This will create a new conda environment called "snakepipes" into which snakePipes is installed. You will now need to create the conda environments needed by the various workflows.

First activate the snakepipes environment with:

``conda activate snakepipes``

Then run the build script with:

``sh build.sh``

You now need to create the various environments required for the pipeline by running:

``snakePipes createEnvs --condaDir ~/.conda/envs/snakepipes``

You do not need to create indices for hg19 or hg38. We provide premade indices stored in a shared location for hg19 and hg38, using the exact fasta and annotation files from the core facility. If you do need additional indices for another organism or genome build, they can be constructed with the following command:

``createIndices --genomeURL <path/URL to your genome fasta> --gtfURL <path/url to genes.gtf> --local -o <output_dir> <name>`` 

Be careful creating indices becuase if you create new indices for hg19 or hg38 , you will change that path in your installation and no longer be using the premade indices. It's best to give any new indices a new name (i.e. hg38_version_x), then they will be stored as a completely different index location the the premade location will remain intact.


Renaming files
-------------
**Note - all of youre sequencing filenames should contain a wz number (i.e. wz3909). Make sure to submit your samples with a wz number in the name or this script will not work.

Before starting a pipeline, it's best to rename your files. The files from the core come with a very long filename (i.e. 5905_25_wz3909_TGACTTCG_S35.bam) and we will shorten this to just the wz number (i.e. wz3909.bam). 

To accomplish this, we have provided an R script above (rename_files.R). This script can either be run from within R, or from the terminal. To run from within R, set your working directory to the folder contaning your files (bam or fastq):

``setwd("/DATA/your.name/your_files/")``

If you copy the script to the same folder as your files, you can run:

``source('rename_files.R')``

Otherwise is can be run from the RStudio script window. 

If you prefer, you can also run from terminal by copying the script into the folder containing your files and running:

``Rscript rename_files.R``

Either way, this will rename all your files and move them into a folder called "rename". All of the files should have been moved into this folder, so if there are any remaining then something went wrong and you should seek help.

Once your files are renamed, you are now ready to proceed with the arropriate pipeline below.

DNA-mapping
-------------

For the DNA-mapping pipeline, the minimum required command is:

``DNA-mapping -i /INPUT/DIR -o /OUTPUT/DIR --local genome_build`` 

The default mapping program is Bowtie2. To use BWA (recommended for Zwart lab ChIP expriments), supply the path to the location of the bwa_mapping.yaml downloaded with this hub. After the renaming step above, all of you files should be in a folder called rename. For example, to run DNA mapping with BWA to hg19, run the following command:

``DNA-mapping -i /PATH/TO/FASTQ/rename -o /PATH/TO/OUTPUT/DIRECTORY --configfile /PATH/TO/snakepipes/bwa_mapping.yaml --local -j 10 --mapq 20 --trim --trim_prg cutadapt --fastqc hg19``

Here, -i specifies the input folder contaning the .fastq.gz files, -o is the output directory, --local runs on the local server and not on a cluster, -j specifies the number of threads, --trim tells the pipeline to trim the reads, --trim_prg tells the pipeline the program used to trim the reads, --fastqc tell it to run fastqc analysis, and finally hg19 specifies the genome build.

ChIP-seq from DNA-mapping pipeline
----------------------------------

The ChIP-seq pipline is designed to take the ouput directly from the DNA-mapping pipeline. The only additional file you will need is a sample_config.yaml file, telling the program your sample names, the control for each sample, and whether to look for broad peaks (i.e. histone marks) or narrow peaks (i.e. transcription factors). See the example sample_config.yaml file above.

If you have run the DNA-mapping pipeline first, then simply run:

``ChIP-seq -d /PATH/TO/DNA-mapping/OUTPUT --local -j 10 --single-end hg19 sample_config.yaml``

Here -d should be the directory with the output of the DNA-mapping pipeline, and it will also direct the output of the ChIP-seq pipeline there. If your samples are not single end then remove the --single-end flag. Also modify the genome_build (i.e. hg19) to suit your purposes).


ChIP-seq from bam files
-----------------------

If you have not run the DNA-mapping pipeline first, then you can still run the pipeline directly from BAM files. In this case,  all of your .bam files shold be renamed in a folder called "rename". You will also need to supply the path to the from_bam.yaml in the snakepipes folder downloaded from this hub. Then run:

``ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bam/rename --configfile /PATH/TO/snakepipes/from_bam.yaml --local -j 10 --single-end hg19 sample_config.yaml``

There will be various folder outputs, including some QC, and the peak files will be in the MACS2 folder. For narrow peaks, the macs2 output will end in ".narrowPeaks", and we have added chr to the chromosome numbers in the file ending in ".chr.narrowPeaks"

Also, creating indices will take some time so you may want to run it in screen to avoid interruptions. (i.e. just add screen -dm before your command, like this: 

``screen -dm ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bam/rename --configfile /PATH/TO/snakepipes/from_bam.yaml --local -j 10 --single-end hg19 sample_config.yaml``

It will look like nothing is happening, but it is running in detached mode and will not be interrupted if your session disconnects. You can see what screens you have running with:

``screen -ls``

If you run screen -ls immediately after executing your screen -dm ChIP-seq... command and you do not see an output for your running screen, then something was wrong with your command (or your environment isn't activated).

Additional Pipelines
-----------------------
The other modules have remained untouched and should work according to the original pipeline `https://github.com/maxplanck-ie/snakepipes`


Finishing up
-------------

When you are finished you should deactivate your conda session to leave the environment with:

``conda deactivate``

This is a good practice so that you don't unintentially alter the environment.

Every time you want to run more analysis you can simply activate the conda environment again with:

``conda activate snakepipes``

All the previously created environments and indices will still be there and you can proceed directly to the pipeline.


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
