FROM ubuntu
MAINTAINER Evan Lezar <evanlezar@gmail.com>

RUN apt-get update && \
	apt-get install -y build-essential

RUN apt-get install -y gfortran \
	libatlas-dev \
	liblapack-dev \
	libarpack2-dev
	
RUN apt-get install -y \
	python \
	python-numpy


COPY src/* /opt/gpu-arpack/src/
COPY test/* /opt/gpu-arpack/test/

WORKDIR /opt/gpu-arpack/src

CMD bash