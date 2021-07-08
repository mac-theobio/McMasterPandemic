FROM r-base:latest

RUN mkdir /home/McMasterPandemic

RUN apt-get update && \
    apt install -y git && \
    apt install -y libcurl4-openssl-dev && \
    apt install -y libssl-dev && \
    apt install -y libxml2-dev && \
    apt install -y pandoc

RUN cd /home && \
    git clone https://github.com/mac-theobio/McMasterPandemic.git

WORKDIR /home/McMasterPandemic

RUN make dependencies

RUN echo "PATH=\"$HOME/bin:$PATH\"" >> /root/.bashrc && \
    . /root/.bashrc && \
    make build-package && \
    make install
