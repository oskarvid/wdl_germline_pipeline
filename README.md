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
As long as you set WORKINGDIR to the path where your wdl, json and conf files are, you don't need to change the settings for their location.  
The .json file assumes that the files are named as follows: GenomeAnalysisTK.jar, picard.jar and cromwell-25.jar.  

You are now ready to run the script:  
```
sh startFromDocker.sh
```
Be adviced that while using the provided toy fastq files the pipeline will crash when it reaches "VariantRecalibration" due to too few samples in the test fastq files. As long as you have a complete set of fastq files the pipeline will finish successfully though.  

### Optimizing  
You should edit the wdl file and change the -t and -nt settings to a value suitable for the machine you'll be running it on. Make sure to change the -Xmx RAM value to something suitable as well. Anything less than 8GB RAM is most likely a bad idea for any full sized dataset.  
Since BaseRecalibrator, PrintReads and HaplotypeCaller are parallelized with scatter gather parallelization, you should edit the application.conf file and change the "concurrent-job-limit" to a value equal to, or less than, the number of available threads your machine has. 

It's a good idea to have at least as many interval lists as you have concurrent jobs running. The current setup has 10 interval lists and 4 concurrent jobs, and as an example our production machine is set to run 18 concurrent jobs with 50 interval lists. Take note that if you have e.g 10 concurrent jobs and 64GB RAM, you need to make sure that each job doesn't consume more than 64/10 = 6.4GB RAM, in a case like that I'd set the -Xmx value to 5 so there'd be some room to spare for the OS. 