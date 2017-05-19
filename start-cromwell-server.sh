# Data file paths
REFERENCE="/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer"
INPUT="/home/oskar/01-workspace/00-temp/wdl_pipeline/data"

# Interval list files
INTERVALS="/home/oskar/01-workspace/00-temp/wdl_pipeline/intervals"

# Tools folder
TOOLS="/home/oskar/01-workspace/00-temp/wdl_pipeline/tools"

# Working directory
WORKINGDIR="/home/oskar/01-workspace/00-temp/wdl_pipeline"

docker run --rm -t -p 8000:8000 \
-v /home/oskar/01-workspace/00-temp/wdl_pipeline/application.conf:/cromwell/application.conf \
-v=$INTERVALS:/intervals \
-v=$WORKINGDIR:/wdl_pipeline \
-v=$TOOLS:/tools \
-v=$REFERENCE:/references \
-v=$INPUT:/data \
oskarv/wdl \
java -jar -Dconfig.file=/cromwell/application.conf \
/tools/cromwell-27.jar server
