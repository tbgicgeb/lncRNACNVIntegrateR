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
    tk-dev

# Install devtools and BiocManager
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages(c('readr', 'dplyr', 'ggpubr', 'survminer', 'PredictABEL', 'enrichR', 'here'))"
RUN R -e "install.packages('BiocManager')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('DESeq2', 'SummarizedExperiment', 'RCircos', 'RegParallel', 'GenomeInfoDb', 'BiocGenerics', 'IRanges', 'GenomicRanges'))"

# Install additional packages from GitHub if necessary
RUN R -e "devtools::install_github('JosephCrispell/basicPlotteR')"

# Copy your package source code into the image
COPY . /usr/local/src/lncRNACNVIntegrateR

# Set the working directory
WORKDIR /usr/local/src/lncRNACNVIntegrateR

# Install your R package
RUN R CMD INSTALL .

# Default command to run R
CMD ["R"]
