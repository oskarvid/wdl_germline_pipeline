# Unofficial development location for a wdl based germline pipeline with scatter gather parallelization

### Disclaimer: The current stable version runs on the b37 decoy reference genome, a version with GRCh38 is in development.

### How to set up and run the pipeline with Docker  
Follow the appropriate instructions on Dockers website for installation instructions: https://docs.docker.com/engine/installation/  
Run this command to download our Docker image for WDL execution:  
```
sudo docker pull oskarv/wdl:latest
```
This docker image has all the dependencies to run our WDL pipeline, including bwa version 0.7.12-r1039. You need to download, GATK-3.7.jar, picard.jar, WDL.jar and cromwell.jar yourself to use this Docker image. Place all of them in the "tools" folder from this repository.  
This article explains how to install GATK and Picard: https://software.broadinstitute.org/gatk/guide/article?id=2899  
And here's instructions for cromwell and wdl: https://software.broadinstitute.org/wdl/userguide/  

You will need to download the reference files too, follow the instructions to access the ftp site [here](https://software.broadinstitute.org/gatk/download/bundle) and download the b37 reference files.

The startFromDocker.sh script assumes that a certain folder structure is in place, and it is recommended that you use the same folder structure unless you want to edit the startFromDocker.sh script and .json file manually to fit your own set up.  
The simplest way to get started is to clone this repository like so:
```
git clone https://github.com/oskarvid/wdl_pipeline.git
```
Now edit the startFromDocker.sh script with your favorite text editor as described below, and then you're good to go!  

Scroll down and change the file paths at "# Data File Paths", "# Tools folder", "# Working directory" and "# Script paths" to point to their respective folders. Do not use a trailing "/" at the end of the file paths, thus "/file/path" is correct wile "/file/path/" is not.  
"REFERENCE" is where your downloaded b37 reference files are.  
"INPUT" is where your input files are, e.g ~/wdl_pipeline/data.  
"TOOLS" is where you placed picard, GATK, cromwell and wdl, e.g ~/wdl_pipeline/tools.  
"WORKINGDIR" is where your output data will get written, e.g ~/wdl_pipeline.  
As long as you set WORKINGDIR to the path where your wdl, json and conf files are, you don't need to change those settings.  
The .json file assumes that the files are named as follows: GenomeAnalysisTK.jar, picard.jar and cromwell-25.jar.  

You are now ready to run the script:  
```
sh startFromDocker.sh
```
Be adviced that while using the provided toy fastq files the pipeline will crash when it reaches "VariantRecalibration" due to too few samples in the test fastq files. As long as you have a complete set of fastq files the pipeline will finish successfully though.  
