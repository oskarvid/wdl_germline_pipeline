echo 'Build Docker Image'
TAG=`date "+%Y%m%d"`
sudo docker build -t oskarv/wdl:$TAG --rm=true .
