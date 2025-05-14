FROM vanallenlab/miniconda:3.12

LABEL description="CAGEcleaner - reducing genomic redundancy in genome mining hit sets"
LABEL org.opencontainers.image.authors="lucas.devrieze@vib.be"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# install an additional perl package (Time::HiRes)
RUN apt-get install -y libmoosex-role-timer-perl unzip

# Install CAGEcleaner and add environment to PATH
RUN conda install -c bioconda -c conda-forge cagecleaner -y
ENV PATH=/usr/miniconda3/envs/cagecleaner/bin/:$PATH
