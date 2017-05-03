FROM ubuntu:16.04
MAINTAINER Oskar Vidarsson <oskar.vidarsson@uib.no>

# Placeholder in case I want to add the tools to the image
#ADD wdltool-0.4.jar /tools/
#ADD cromwell-0.21.jar /tools/
#ADD GenomeAnalysisTK.jar /tools/GenomeAnalysisTK-3.6/

# Install wget
RUN apt-get update && apt-get install -y \
bwa \
wget \
&& rm -rf /var/lib/apt/lists/*

# Clean
RUN rm -rf /usr/share/locale/ /usr/share/man/ /root/.cache

# Install Java 8
RUN mkdir /opt/jdk && cd /opt && \
wget --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u5-b13/jdk-8u5-linux-x64.tar.gz && \
tar -zxf jdk-8u5-linux-x64.tar.gz -C /opt/jdk && \
rm jdk-8u5-linux-x64.tar.gz && \
update-alternatives --install /usr/bin/java java /opt/jdk/jdk1.8.0_05/bin/java 100 && \
update-alternatives --install /usr/bin/javac javac /opt/jdk/jdk1.8.0_05/bin/javac 100


