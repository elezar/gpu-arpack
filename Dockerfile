FROM ubuntu
MAINTAINER evanlezar@gmail.com

COPY src /opt/gpu-arpack/
COPY test /opt/gpu-arpack/

WORKDIR /opt/gpu-arpack/src

RUN make

CMD bash