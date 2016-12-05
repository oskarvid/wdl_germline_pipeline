# Maintainer: Oskar Vidarsson <oskar.vidarsson@uib.no>

#Container ID
IMAGE_ID="f64453d5ba50"

#User ID (Choose one of CUST_USERID definition)
USERID=`id -u`
GROUPID=`id -g`
#CUST_USERID="-u=501:501"
CUST_USERID="-u=$USERID:$GROUPID"

# Data file paths
REFERENCE="/path/to/reference/files"
INPUT="/path/to/wdl_pipeline/data"

# Tools folder
TOOLS="/path/to/wdl_pipeline/tools"

# Working directory
WORKINGDIR="/path/to/wdl_pipeline"

# Script paths
WDLSCRIPT="germlinevarcall.wdl"
WDLJSON="germlinevarcall.json"

# Job
JAVA="java -jar /tools/cromwell-0.21.jar run $WDLSCRIPT $WDLJSON"

# For Debug
#sudo docker run -ti --rm $CUST_USERID -v=/$TOOLS:/tools -v=$REFERENCE:/reference -v=$DATA:/data -w=/data/ $IMAGE_ID bash
# Start without --rm to be able to commit changes of image
#docker run --rm -ti $CUST_USERID -v=$WORKINGDIR:/wdl_pipeline -v=$TOOLS:/tools -v=$REFERENCE:/references -v=$INPUT:/data -w=/wdl_pipeline $IMAGE_ID bash

# For Running
sudo docker run -t --rm $CUST_USERID -v=$WORKINGDIR:/wdl_pipeline -v=$TOOLS:/tools -v=$REFERENCE:/references -v=$INPUT:/data -w=/wdl_pipeline --cpuset-cpus="0-3" -c 512  $IMAGE_ID sh -c "$JAVA"
