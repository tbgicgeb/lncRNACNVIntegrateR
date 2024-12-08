# Start from the Rocker R image
FROM rocker/r-ver:4.1.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libxt-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    tcl8.6-dev \
    tk8.6-dev \
    tk-dev \
    wget

# Install devtools and a specific version of BiocManager
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('BiocManager', version = '1.30.25')"

# Install general R packages with specific versions
RUN R -e "install.packages(c('readr' = '2.1.5', 'dplyr' = '1.1.4', 'ggpubr' = '0.6.0', 'survminer' = '0.4.9', 'PredictABEL' = '1.2-4', 'enrichR' = '3.2', 'here' = '1.0.1'))"

# Set Bioconductor version to 3.13 and install specific packages with versions
RUN R -e "BiocManager::install(c('DESeq2' = '1.32.0', 'SummarizedExperiment' = '1.24.0', 'RCircos' = '1.2.2', 'GenomeInfoDb' = '1.30.1', 'BiocGenerics' = '0.40.0', 'IRanges' = '2.28.0', 'GenomicRanges' = '1.46.1', 'RegParallel' = '1.10.0'))"

# Install additional packages from GitHub if necessary
RUN R -e "devtools::install_github('JosephCrispell/basicPlotteR')"

# Copy your package source code into the image
COPY . /usr/local/src/lncRNACNVIntegrateR

# Set the working directory
WORKDIR /usr/local/src/lncRNACNVIntegrateR

# Install your R package
RUN R CMD INSTALL .

# Default command to run R
CMD ["R", "--no-save"]

