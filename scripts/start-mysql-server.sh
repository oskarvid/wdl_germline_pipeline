sudo docker run --rm -t \
--env="MYSQL_ROOT_PASSWORD=cromwell" \
--env="MYSQL_DATABASE=cromwell_db" \
-p 3308:3306 \
-v `pwd`/cromwell-mysql/mysql/init:/docker-entrypoint-initdb.d \
-v `pwd`/cromwell-mysql/mysql/data:/var/lib/mysql mysql:5.7
