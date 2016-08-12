FROM debian:jessie
MAINTAINER Jake Hagen <jake.hagen@opmbx.org>
RUN echo deb http://cloud.r-project.org/bin/linux/debian jessie-cran3/ >> /etc/apt/sources.list
RUN apt-get update && apt-get install -y --force-yes \
		r-base=3.3.1* \
		r-base-dev=3.3.1* \
		r-base-dev=3.3.1* \
		r-cran-xml \
		libcurl4-openssl-dev
COPY install_script.R /tmp/
COPY diffExpression.R /usr/local/bin/diffExpression.R
RUN Rscript /tmp/install_script.R
ENTRYPOINT ["Rscript", "/usr/local/bin/diffExpression.R"]
