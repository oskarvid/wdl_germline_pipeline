curl -X POST --header "Accept: application/json" "http://localhost:8001/api/workflows/v1/batch" \
	-F workflowSource=@germlinevarcall_bwa_loop.wdl \
	-F workflowInputs=@germlinevarcall-hg38.json \
	-F workflowOptions=@cromwell-mysql/cromwell/workflow-options/workflowoptions.json
