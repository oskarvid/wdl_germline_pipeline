# WDL based best practices germline variant calling pipeline with scatter gather parallelization

### How to set up and run the pipeline with Docker  
Follow the appropriate instructions on Dockers' website for installation instructions: https://docs.docker.com/engine/installation/  
Run this command to download the Docker image for WDL execution:  
```
sudo docker pull oskarv/wdl
```
And download the MySQL docker image too:
```
sudo docker pull mysql:5.7
```
The oskarv/wdl docker image has all the software dependencies to run the WDL pipeline, including bwa version 0.7.12-r1039 and picard 2.9.1. You need to download GATK-3.7 and cromwell-27.jar yourself to use this Docker image. At the time of writing this, cromwell 27 isn't available for regular download, but can be either buil manually from the [cromwell repository](https://github.com/broadinstitute/cromwell) or copied from a downloaded broadinstitute/cromwell:develop docker image. Cromwell 26 or any other earlier version will not work with the supplied version of the application.conf file. Once you've downloaded GATK and cromwell, place them in the "tools" folder from this repository.  
This article explains how to download GATK: https://software.broadinstitute.org/gatk/guide/article?id=2899  
And here's instructions for wdl: https://software.broadinstitute.org/wdl/userguide/  
N.B: You don't strictly need wdl.jar to run the pipeline in this repository, but it's usefull for creating json files and such.


The current version of the pipeline runs on the b37 reference genome, but it will be upgraded to GRCh38 in the future. You need to download the b37 reference files from here: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/, if it asks you for a password, just hit enter, there is no password. Or you can use your preferred reference version, and in that case you need to edit the json files and wdl files.  

### Configuring the scripts and running the pipeline
Let's get started with setting up the scripts so we can run the pipeline.  
Fire up a terminal window and clone this repository:
```
git clone https://github.com/oskarvid/wdl_pipeline.git
```
As long as you run the pipeline from this folder, the only thing you need to edit is the start-cromwell.sh script and edit the directory path for the reference file location.   

But if you want to use e.g another input file location, open the start-cromwell.sh script and edit the "INPUT=" line with the correct path to your input data.  

You are now ready to start the MySQL server:  
```
sh start-mysql-server.sh
```
Once it's up and running you can start the cromwell server:
```
sh start-cromwell-server.sh
```
And when cromwell is up and running, you can run the pipeline with the test data like so:
```
sh start-pipeline.sh
```
Be adviced that when you use the provided toy fastq files and run the germlinevarcall.wdl pipeline, it will crash when it reaches "VariantRecalibration" due to too few reads in the test fastq files. As long as you have a complete set of fastq files the pipeline will finish successfully.  

If you want to kill the pipeline, you can use the kill-wdl-job.sh script, or go to localhost:8000 and use the control panel there.  
When you've done some testing you can use the clean.sh script to remove the MySQL database, cromwell-executions folder as well as the cromwell-workflow-logs folder.  

### But what do the start scripts do?

Without going into too much detail, here's a description of what the start scripts do.  
1. start-mysql-server.sh takes a .sql file as input which creates a database called cromwell_db, as well as a user named cromwell, and gives the user rights to read and write to cromwell_db. This init script is located in cromwell-mysql/mysql/init. It's important that you always start the MySQL database first, otherwise cromwell won't start.
2. start-cromwell.sh mounts a set of directories as volumes inside the docker container, starts the cromwell server and connects to the MySQL database. Cromwell also starts a web server that can be reached by going to localhost:8000, from there you have access to a basic GUI that lets you start jobs, stop jobs, view metadata and more.  
3. start-pipeline.sh sends a curl command with the wdl script, input json and workflow options json file to the cromwell api at localhost:8000. This tells cromwell to start the pipeline and do its thing. The output files are put in the folder "cromwell-executions". 