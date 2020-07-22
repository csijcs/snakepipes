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

Begin by logging into harris server and entering the terminal environment.

First initialize conda in your environment with:

``conda init``

Ensure conda is properly installed by running:

``conda --version``

Find the path to the conda your are using with:

``which conda``

For harris server, your conda path should be:

``/opt/anaconda/bin/conda``

Once you have conda working and the proper path, configure the directory for pkgs to be installed with:

``/opt/anaconda/bin/conda config --add pkgs_dirs ~/.conda/pkgs/``

Then clone this repository into your desired location with:

``git clone https://github.com/csijcs/snakepipes.git``

Change directory into the snakepipes folder with:

``cd snakepipes``

First, create a new conda environment called "snakepipes" into which snakePipes is installed with:

``/opt/anaconda/bin/conda env create --file snakepipes.yaml``

Next, activate the snakepipes environment with:

``conda activate snakepipes``

Then run the build script with:

``sh build.sh``

You now need to create the various environments required for the pipelines by running:

``snakePipes createEnvs --condaDir ~/.conda/envs/snakepipes``

You do not need to create indices for hg19 or hg38. We provide premade indices stored in a shared location for hg19 and hg38, using the exact fasta and annotation files from the core facility. If you do need additional indices for another organism or genome build, seek assistance from a computational expert.

Renaming files
-------------
**Note - all of youre sequencing filenames should contain a wz number (i.e. wz3909). Make sure to submit your samples with a wz number in the name or this script will not work. If there are two samples with the same wz number (i.e. same sample split across two lanes) with second file will be renamed wzNUMBER_2 (i.e. wz3909_2). If there are more than two (not likely, but possible), it will give an error and not rename your additional files. If you do actually have more than two (i.e. split across more than two lanes), seek professional help.

Before starting a pipeline, it's best to rename your files. The files from the core come with a very long filename (i.e. 5905_25_wz3909_TGACTTCG_S35.bam) and we will shorten this to just the wz number (i.e. wz3909.bam). 

To accomplish this, we have provided an R script above (rename_files.R). This script can either be run from within R, or from the terminal. To run from within R, set your working directory to the folder contaning your files (bam or fastq):

``setwd("/DATA/your.name/your_files/")``

If you copy the script to the same folder as your files, you can run:

``source('rename_files.R')``

Otherwise is can be run from the RStudio script window. 

If you prefer, you can also run from terminal by copying the script into the folder containing your files and running:

``Rscript rename_files.R``

Either way, this will rename all your files and move them into a folder called "rename". All of the files should have been moved into this folder, so if there are any remaining then something went wrong and you should seek help.

Once your files are renamed, you are now ready to proceed with the appropriate pipeline below.

DNA-mapping
-------------

For DNA mapping, we generally recommend using BWA. To do this, supply the path to the location of the bwa_mapping.yaml downloaded with this hub. After the renaming step above, all of your fastq files should be in a folder called rename. Be sure you know the appropriate genome build for your project (i.e. hg19 or hg38). For example, to run DNA mapping with BWA to hg19, run the following command:

``DNA-mapping -i /PATH/TO/FASTQ/rename -o /PATH/TO/OUTPUT/DIRECTORY --configfile /PATH/TO/snakepipes/bwa_mapping.yaml --local -j 10 --mapq 20 --trim --trim_prg cutadapt --fastqc hg19``

Here, -i specifies the input folder contaning the fastq files, -o is the output directory of your choosing, and  hg19 specifies the genome build (adjust to hg38 as appropriate for your project). The rest of the parameters should not be altered for standard ChIP-seq experiments.

**Note - Previous projects as well as many existing projects in the Zwart lab have been mapped using the bwa-backtrack algorithm. For legacy reasons, if you need your peakcalling results to match EXACTLY to previous results, we recommend using the bam files supplied by the core and taking them through the ChIP-seq from bam pipeline below. The BWA option for this DNA-mapping pipeline uses the bwa-mem algorithm, which will produce very similar but not exactly the same results.  

If, for purposes other than Zwart lab ChIP expirements, you would like to map with Bowtie, simply remove the --configfile /PATH/TO/snakepipes/bwa_mapping.yaml from the command.


ChIP-seq from DNA-mapping pipeline
----------------------------------

The ChIP-seq pipline is designed to take the ouput directly from the DNA-mapping pipeline. The only additional file you will need is a sample_config.yaml file, telling the program your sample names, the control for each sample, and whether to look for broad peaks (i.e. histone marks) or narrow peaks (i.e. transcription factors). See the example sample_config.yaml file above.

If you have run the DNA-mapping pipeline first, then simply run:

``ChIP-seq -d /PATH/TO/DNA-mapping/OUTPUT --local -j 10 --single-end hg19 sample_config.yaml``

Here -d is the directory with the output of the DNA-mapping pipeline, and it will also direct the output of the ChIP-seq pipeline there. 

**Note - The new projects should be getting mapped to the hg38 genome build, while ongoing projects that were previously mapped to hg19 should stay with hg19. Ensure you are not mixing hg38 and hg19 in your project or the results will be inconsistent.  

**Note - Most, if not all, Zwart lab ChIP experiments will be single-end. If you have paired-end reads from a collaborator or publically available dataset, you will need to supply the paired_end_from_bam.yaml file instead, and remove the --single-end option.

ChIP-seq from bam files
-----------------------

If you have not run the DNA-mapping pipeline, you can run the ChIP-seq pipeline directly from BAM files. In this case, all of your .bam files should be renamed in a folder called "rename". You will also need to supply the path to the from_bam.yaml in the snakepipes folder downloaded wit this hub. For single-end reads the command to run is:

``ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bam/rename --configfile /PATH/TO/snakepipes/from_bam.yaml --local -j 10 --single-end hg19 sample_config.yaml``

There will be various folder outputs, including some QC, and the peak files will be in the MACS2 folder. For narrow peaks, the macs2 output will end in ".narrowPeaks", and we have added chr to the chromosome numbers in the file ending in ".chr.narrowPeaks" for your convenience.

**Note - The new projects should be getting mapped to the hg38 genome build, while ongoing projects that were previously mapped to hg19 should stay with hg19. Ensure you are not mixing hg38 and hg19 in your project or the results will be inconsistent.  

**Note - Most, if not all, Zwart lab ChIP experiments will be single-end. If you have paired-end reads from a collaborator or publically available dataset, you will need to supply the paired_end_from_bam.yaml file instead, and remove the --single-end option.

Running Pipelines in screen
----------------------------
Running pipelines will take some time, so you will want to run in screen to avoid interruptions. To do this, just add screen -dm before your command, like this: 

``screen -dm ChIP-seq -d /PATH/TO/OUTPUT/DIR --fromBam /PATH/TO/bam/rename --configfile /PATH/TO/snakepipes/from_bam.yaml --local -j 10 --single-end hg19 sample_config.yaml``

It will look like nothing is happening, but it is running in detached mode and will not be interrupted if your session disconnects. Furthermore, it will disconnect automatically when it is finished. You can see what screens you have running with:

``screen -ls``

If you run screen -ls immediately after executing your screen -dm ChIP-seq... command and you do not see an output for your running screen, then something went wrong (or your environment isn't activated).

Additional Pipelines
-----------------------
The other modules have remained untouched and should work according to the original pipeline `https://github.com/maxplanck-ie/snakepipes`


Finishing up
-------------

When you are finished you should deactivate your conda session to leave the environment with:

``conda deactivate``

This is a good practice so that you don't unintentially alter the environment. 

Never install anything else within your snakepipes environment.

Every time you want to run more analysis you can simply activate your environment again with:

``conda activate snakepipes``

All the previously created environments and indices will still be there and you can proceed directly to the pipelines.


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
