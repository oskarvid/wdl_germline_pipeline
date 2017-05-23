curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
	-F wdlSource=@germlinevarcall.wdl \
	-F workflowInputs=@germlinevarcall.json \
	-F workflowOptions=@cromwell-mysql/cromwell/workflow-options/workflowoptions.json