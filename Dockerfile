FROM ubuntu:22.04
COPY . /app
RUN apt-get update \
    && apt-get install -y \
        git g++ cmake \
    && apt-get clean

RUN apt-get install -y libboost-all-dev \
        libtiff-dev \
        libgeotiff-dev \
        install pdal \
    && rm -rf /var/lib/apt/lists/*

RUN cd /app

WORKDIR /app

CMD ["bash"]