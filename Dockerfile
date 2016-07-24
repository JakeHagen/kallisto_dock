FROM debian:jessie
MAINTAINER Jake Hagen <jake.hagen@opmbx.org>
RUN apt-get update && apt-get install -y \
		build-essential \
		cmake \
		git \
		libhdf5-dev \
		zlib1g-dev
WORKDIR /docker
RUN git clone https://github.com/pachterlab/kallisto.git
WORKDIR /docker/kallisto
RUN mkdir build
WORKDIR /docker/kallisto/build
RUN cmake .. && \
	make && \
	make install
Volume ["/home/jake/rna_data/"]
ENTRYPOINT ["kallisto"]
