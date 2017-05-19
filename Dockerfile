FROM ubuntu:16.04
MAINTAINER Oskar Vidarsson <oskar.vidarsson@uib.no>

# Install wget and bwa
RUN apt-get update && apt-get install -y \
bwa \
wget && \
rm -rf /var/lib/apt/lists/* && \
rm -rf /usr/share/locale/ /usr/share/man/ /root/.cache

# Download picard
RUN mkdir /Jar && \
wget https://github.com/broadinstitute/picard/releases/download/2.9.1/picard.jar -O /Jar/picard.jar

# Install Java 8
RUN mkdir /opt/jdk && cd /opt && \
wget --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u5-b13/jdk-8u5-linux-x64.tar.gz && \
tar -zxf jdk-8u5-linux-x64.tar.gz -C /opt/jdk && \
rm jdk-8u5-linux-x64.tar.gz && \
update-alternatives --install /usr/bin/java java /opt/jdk/jdk1.8.0_05/bin/java 100 && \
update-alternatives --install /usr/bin/javac javac /opt/jdk/jdk1.8.0_05/bin/javac 100
