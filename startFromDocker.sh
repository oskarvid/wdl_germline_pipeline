# Maintainer: Oskar Vidarsson <oskar.vidarsson@uib.no>

#Container ID
IMAGE_ID="oskarv/wdl:latest"

#User ID (Choose one of CUST_USERID definition)
USERID=`id -u`
GROUPID=`id -g`
CUST_USERID="-u=$USERID:$GROUPID"

# Data file paths
REFERENCE="/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer"
INPUT="/home/oskar/01-workspace/00-temp/wdl_pipeline/data"

# Tools folder
TOOLS="/home/oskar/01-workspace/00-temp/wdl_pipeline/tools"

# Working directory
WORKINGDIR="`pwd`"

# Script paths
WDLSCRIPT="germlinevarcall.wdl"
WDLJSON="germlinevarcall.json"
CROMWELLSETTINGS="application.conf"
OPTIONS="workflowoptions.json"

# Job
JAVA="java -Dconfig.file=/wdl_pipeline/application.conf -jar \
/tools/cromwell-27.jar run $WDLSCRIPT $WDLJSON"

# For Debug
#docker run --rm -ti $CUST_USERID -v=$WORKINGDIR:/wdl_pipeline -v=$TOOLS:/tools -v=$REFERENCE:/references -v=$INPUT:/data -w=/wdl_pipeline $IMAGE_ID bash

# For Running
sudo docker run -t --rm $CUST_USERID -v=$WORKINGDIR:/wdl_pipeline \
-v=$TOOLS:/tools -v=$REFERENCE:/references -v=$INPUT:/data -w=/wdl_pipeline \
--cpuset-cpus="0-3" -c 512  $IMAGE_ID sh -c "$JAVA"
