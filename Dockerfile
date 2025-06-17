FROM ubuntu:22.04
COPY . /workspace

RUN apt-get update \
    && apt-get install -y git g++ gcc make gdb \
    && apt-get install -y pdal libpdal-dev \
    && apt-get clean

WORKDIR /workspace

CMD ["bash"]