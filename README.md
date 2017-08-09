# Germline Variant Calling Pipeline Based On WDL, Cromwell and MySQL

There are two ways to run this pipeline, either with the ad hoc method that will let you run it without support for restarting stopped pipelines, reusing previous results or proper process limitation for resource control. Or with all bells and whistles that requires running Cromwell as a server as well as running a MySQL database.  
The first instructions will describe how to run the basic and simpler version of the pipeline, keep in mind that it's not recommended for production settings (and the moment it's running GATK 4 beta which is only suitable for research or testing purposes).  

## Setup instructions for the basic pipeline execution mode
To run the simple version of the pipeline you'll use the start-FromDocker.sh script that is located in the start-scripts folder. You need to edit it according to where your reference files, input files etc are. Also edit the germlinevarcall.json file with the correct file names, if you don't change the docker mount points in the start-FromDocker.sh script, you should be fine if you just edit the file names in the .json file.  
Finally, you need to download cromwell with the dl-docker.sh script, and obviously you need to have docker installed. When everything it set up, simply run "sh start-FromDocker.sh"  

## Setup instructions for the full pipeline execution mode
To run the pipeline with all bells and whistles, these three scripts need to be executed in the following order.  
1. sh start-mysql-server.sh  
2. sh start-cromwell-server.sh  
3. sh start-pipeline.sh  

If you simply run the start-mysql-server.sh script it will download the required mysql:5.7 docker image automatically. You might want to edit the username and password, the default is set to "cromwell" for both.  

The second step is setting up the start-cromwell-server.sh script. As long as you run the pipeline from the "wdl_pipeline" directory you only need to edit the REFERENCE, which is the directory where your reference files are located. You also need to run "sudo docker ps -a" and copy the container ID of the MySQL docker container, e.g "3da13d9f19b0", then run "docker inspect 3da13d9f19b0 | grep IPA" and copy the IP address and paste it in the application.conf file that is located at cromwell-mysql/cromwell/config/. Towards the bottom of the file there is an IP address the points to the MySQL database, replace it with your copied IP address.  
Before you can run the start-cromwell.sh script, you need to run the dl-cromwell.sh script to download cromwell. It will automatically get downloaded to the tools directory, and now you can finally run "sh start-scripts/start-cromwell-server.sh".

Now that the MySQL and cromwell servers are up and running you can run the start-pipeline.sh script. It will use the fastq files that are in the data directory. Keep in mind that the pipeline will fail when it reaches the VariantRecalibration step since there isn't enough variants, complete fastq files obviously work though.

## The finer details

The pipeline is compartmentalised into five major blocks. The MySQL server, the Cromwell server, the WDL script, the json file with the reference file- and a few other inputs (not to be confused with the json file with the workflow options, more on that later), and the tsv file which defines the attributes of the fastq files. We've already looked closer at the MySQL- and Cromwell server scripts, so we'll continue with the WDL script.  
The WDL script contains the pipeline instructions, this is where you go to change how the pipeline is executed.  
The json file with the file paths to the reference files, interval files, string definitions etc, is where you change the input files and some of the input strings for the WDL script.  
The template_sample_manifest.tsv file in the intervals directory is where you put the file paths to your fastq files, as well as where you add information about sample names, flowcells, lanes, libraries and platform. This information is needed for correct read group creation for the bam files. BaseRecalibrator and MarkDuplicates benefit from having lane- and library information. If you put your input files in the wdl_pipeline/data directory you don't need to edit the file path, just edit the name from "test_R1.fastq.gz" to whatever your file name is.  

There are two more configuration files of interest, namely the workflowoptions.json file and the application.conf files.  
The workflowoptions.json file is used to set parameters for call caching, and there's a place holder for the final workflow output directory, but it doesn't work as intended at the time of writing this. It can be expanded with more options.  
The final config file of interest is the application.conf file in the cromwell-mysql/cromwell/config/ directory. This is the file that is used for setting options for cromwell. You can edit the number of concurrently running jobs, i.e affect how many processes are created during scatter gather execution. It is set to 2 by default only for testing purposes, increase the number to something suitable for your server or laptop. Keep in mind that the tools are very RAM hungry, and each process is standalone, so you need to balance the -Xmx java setting in the WDL script with the number of concurrently running jobs.

There are way more things that I could elaborate on, but this should be enough for most users. A planned feature is support for the GenomeDB database to take advantage of its speed enhancements. None of the tools are parallelized except for bwa, but the ones that will support it once GATK 4 goes stable will get it, most likely GenotypeGVCFs, VariantRecalibration and ApplyVQSR.
