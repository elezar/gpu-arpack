FROM ubuntu
MAINTAINER Evan Lezar <evanlezar@gmail.com>

RUN apt-get update && \
	apt-get install -y \
		build-essential \
		gfortran \
		libatlas-dev \
		liblapack-dev \
		libarpack2-dev \
		python \
		python-numpy && \
	apt-get clean &&\
	rm -f /etc/apt/source.list.d/*


COPY src/* /opt/gpu-arpack/src/
COPY test/* /opt/gpu-arpack/test/

WORKDIR /opt/gpu-arpack/src

RUN make

WORKDIR /opt/gpu-arpack/test

ENTRYPOINT ["python", "python_driver.py"]
