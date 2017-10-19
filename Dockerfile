FROM debian:jessie

MAINTAINER jianglikun@ds.quantibio.com

RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    wget \
    make \
    less \
    build-essential \
    python-dev \
    ca-certificates
RUN wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py && pip install --upgrade ete3 && easy_install -U ete3 && pip install six
ADD ./script /
ADD assembly.id.txt /
ADD 4904MG.info.txt /
