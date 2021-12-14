FROM r-base:latest

RUN mkdir /home/McMasterPandemic

RUN apt-get update && \
    # git, libcurl4-openssl-dev, libssl-dev **not required** on https://hub.docker.com/r/rocker/rstudio
    apt install -y git && \
    apt install -y libcurl4-openssl-dev && \
    apt install -y libssl-dev && \
    # libxml2-dev, pandoc **required** on https://hub.docker.com/r/rocker/rstudio
    apt install -y libxml2-dev && \
    apt install -y pandoc

RUN cd /home && \
    git clone https://github.com/mac-theobio/McMasterPandemic.git

RUN echo "force layer construction"

WORKDIR /home/McMasterPandemic

RUN make dependencies

RUN git fetch origin tmb

RUN git checkout tmb

RUN echo "PATH=\"$HOME/bin:$PATH\"" >> /root/.bashrc && \
    . /root/.bashrc && \
    make install BUILDARGS="--no-build-vignettes "
