FROM ubuntu:16.04
MAINTAINER Oskar Vidarsson <oskar.vidarsson@uib.no>

# Install wget and bwa
RUN apt-get update && apt-get install -y \
bwa \
wget \
unzip \
python-minimal && \
rm -rf /var/lib/apt/lists/* && \
rm -rf /usr/share/locale/ /usr/share/man/ /root/.cache

# Download picard
RUN mkdir /Jar && \
wget https://github.com/broadinstitute/gatk/releases/download/4.beta.3/gatk-4.beta.3.zip -O /Jar/gatk4.zip && \
unzip -q /Jar/gatk4.zip -d /Jar/ && \
mv /Jar/gatk*/gatk* /Jar/ && \
rm -r /Jar/gatk*/ /Jar/gatk4.zip

# Install Java 8
RUN mkdir /opt/jdk && cd /opt && \
wget --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz && \
tar -zxf jdk-8u131-linux-x64.tar.gz -C /opt/jdk && \
rm jdk-8u131-linux-x64.tar.gz && \
update-alternatives --install /usr/bin/java java /opt/jdk/jdk1.8.0_131/bin/java 100 && \
update-alternatives --install /usr/bin/javac javac /opt/jdk/jdk1.8.0_131/bin/javac 100
