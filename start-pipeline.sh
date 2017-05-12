curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
	-F wdlSource=@germlinevarcall-testing-bwa-scatter.wdl \
	-F workflowInputs=@germlinevarcall-server-bwa-scatter.json \
	-F workflowOptions=@workflowoptions.json