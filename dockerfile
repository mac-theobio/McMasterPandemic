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

WORKDIR /home/McMasterPandemic

RUN echo "force layer construction"

RUN make dependencies

RUN echo "PATH=\"$HOME/bin:$PATH\"" >> /root/.bashrc && \
    . /root/.bashrc && \
    make build-package && \
    make install
