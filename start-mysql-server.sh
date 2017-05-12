sudo docker run --rm -t \
--env="MYSQL_ROOT_PASSWORD=cromwell" \
--env="MYSQL_DATABASE=cromwell_db" \
-p 3307:3306 \
-v /home/oskar/01-workspace/00-temp/wdl_pipeline/cromwell-mysql/mysql/init:/docker-entrypoint-initdb.d \
-v /home/oskar/01-workspace/00-temp/wdl_pipeline/cromwell-mysql/mysql/data:/var/lib/mysql mysql:5.7