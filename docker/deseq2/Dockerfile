FROM bioconductor/release_base

MAINTAINER John Vivian, jtvivian@gmail.com

RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('BiocParallel',ask=FALSE)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('DESeq2',ask=FALSE)"

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["Rscript"]
