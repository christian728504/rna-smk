FROM python:3.12-slim-bookworm

ENV PYTHONUNBUFFERED=1
ENV HOME=/root
ENV CXXFLAGS="-O3"
ENV CFLAGS="-O3"

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    wget \
    uuid-runtime \
    lz4 \
    cmake \
    bedtools \
    rsync \
    libhts3 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    pkg-config \
    autoconf \
    lbzip2 \
    libexpat1-dev \
    libncurses5-dev \
    curl \
    libfreetype6-dev \
    libfontconfig1-dev

RUN apt-get update && apt-get install -y locales && \
    sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

RUN wget http://opengene.org/fastp/fastp && \
    mv fastp /usr/local/bin/ && \
    chmod a+x /usr/local/bin/fastp

RUN wget -O /tmp/pigz.tar.gz https://zlib.net/pigz/pigz.tar.gz && \
    cd /tmp && \
    tar -xzf pigz.tar.gz && \
    cd pigz && \
    make && \
    cp pigz /usr/local/bin/ && \
    cp unpigz /usr/local/bin/ && \
    rm -rf /tmp/pigz*

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
RUN echo 'source /root/.cargo/env' >> /root/.bashrc
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo install bigtools
RUN cp /root/.cargo/bin/bigtools /usr/local/bin/

RUN wget https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 -O /opt/samtools-1.22.tar.bz2 && \
    tar -xvf /opt/samtools-1.22.tar.bz2 -C /opt && \
    cd /opt/samtools-1.22 && \
    ./configure && \
    make && \
    make install

RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 -O /opt/htslib-1.22.tar.bz2 && \
    tar -xvf /opt/htslib-1.22.tar.bz2 -C /opt && \
    cd /opt/htslib-1.22 && \
    ./configure && \
    make && \
    make install

RUN curl -L https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz | tar -xzf - -C /opt && \
    ln -s /opt/STAR-2.7.11b/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR-2.7.11b/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong
    
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz | tar -xzf - -C /usr/local/bin

RUN git clone https://github.com/deweylab/RSEM.git /opt/RSEM && \
    cd /opt/RSEM && \
    make && \
    make install

RUN git clone https://github.com/DaehwanKimLab/hisat2.git /opt/hisat-3n && \
    cd /opt/hisat-3n && \
    git checkout -b hisat-3n origin/hisat-3n && \
    make
ENV PATH="/opt/hisat-3n:${PATH}"

RUN git clone https://github.com/christian728504/mudskipper.git /opt/mudskipper && \
    cd /opt/mudskipper && \
    cargo build --release && \
    cp target/release/mudskipper /usr/local/bin/

RUN apt update && apt install -y python3 python3-pip

RUN pip install --no-cache-dir \
    "beautifulsoup4>=4.13.4" \
    "ipykernel>=6.29.5" \
    "jinja2>=3.1.6" \
    "lxml>=5.4.0" \
    "multiprocess>=0.70.18" \
    "pandas>=2.3.0" \
    "polars>=1.30.0" \
    "pysam>=0.23.3" \
    "qc-utils>=20.9.1" \
    multiqc
    
CMD ["/bin/bash"]
