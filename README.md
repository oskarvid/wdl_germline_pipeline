# WDL based best practices germline variant calling pipeline with scatter gather parallelization

# These instructions are not up to date, if you don't know what you're doing you're gonna have a bad time, you have been warned. But try doing something with the start-*.sh scripts, use some other wdl, json, workflowoptions and applications.conf files and you should be good to go though ;)  

### How to set up and run the pipeline with Docker  
Follow the appropriate instructions on Dockers website for installation instructions: https://docs.docker.com/engine/installation/  
Run this command to download our Docker image for WDL execution:  
```
sudo docker pull oskarv/wdl
```
This docker image has all the dependencies to run our WDL pipeline, including bwa version 0.7.12-r1039. You need to download, GATK-3.7 and cromwell-27.jar yourself to use this Docker image. Place all of them in the "tools" folder from this repository.  
This article explains how to install GATK and Picard: https://software.broadinstitute.org/gatk/guide/article?id=2899  
And here's instructions for cromwell and wdl: https://software.broadinstitute.org/wdl/userguide/ (instructions unclear for cromwell-27.jar at the moment, can be copied from the broadinstitute/cromwell docker image)  

  
The simplest way to get started is to clone this repository like so:
```
git clone https://github.com/oskarvid/wdl_pipeline.git
```
Now edit the startFromDocker.sh script as described below, and then you're good to go!  

First of all, run this:  
```
sudo docker images
```
and copy the image ID of your "oskarv/wdl" image, open startFromDocker.sh with a text editor and paste the image id into "IMAGEID=" at the top of the file.  

Scroll down and change the file paths at "# Data File Paths", "# Tools folder", "# Working directory" and "# Script paths" to point to their respective folders. Do not use a trailing "/" at the end of the file paths, thus "/file/path" is correct wile "/file/path/" is not.  
"REFERENCE" is where your reference files are.  
"DATA" is where your input files are, e.g ~/wdl_pipeline/data.  
"WORKINGDIR" is where your output data will get written, e.g ~/wdl_pipeline.  
"TOOLS" is where you placed picard, GATK, cromwell and wdl, e.g ~/wdl_pipeline/tools.  
The .json file assumes that the files are named as follows: GenomeAnalysisTK.jar (i.e GATK 3.6), GATK-4.jar, picard.jar, cromwell-0.21.jar and wdltool-0.4.jar.  

You are now ready to run the script:  
```
sh startFromDocker.sh
```
Be adviced that while using the provided toy fastq files the pipeline will crash when it reaches "VariantRecalibration" due to too few samples in the test fastq files. As long as you have a complete set of fastq files the pipeline will finish successfully though.  
