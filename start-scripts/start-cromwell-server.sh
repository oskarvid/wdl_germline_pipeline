# Data file paths
REFERENCE="/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer"
INPUT="/home/oskar/"
#04-pipelines/GATK-Ghislain/GermlineVarCall-Workflow/Samples/"

# Interval list files
INTERVALS="`pwd`/intervals"

# Tools folder
TOOLS="`pwd`/tools"

# Application config
APPCONFIG="`pwd`/cromwell-mysql/cromwell/app-config/"

# Working directory
WORKINGDIR="`pwd`"

docker run --rm -t -p 8000:8000 \
-v=`pwd`/cromwell-executions:/cromwell-executions \
-v=`pwd`/cromwell-workflow-logs:/cromwell-workflow-logs \
-v=$APPCONFIG:/cromwell \
-v=$INTERVALS:/intervals \
-v=$WORKINGDIR:/wdl_pipeline \
-v=$TOOLS:/tools \
-v=$REFERENCE:/references \
-v=$INPUT:/data \
oskarv/wdl:gatk4 \
java -jar -Dconfig.file=/cromwell/application.conf \
/tools/cromwell-28_2.jar server
