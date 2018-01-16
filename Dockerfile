FROM sjackman/linuxbrew-core
MAINTAINER Wen-Wei Liang <wenwiliang@gmail.com>

USER root
RUN apt-get update
RUN apt-get install -y pkg-config autoconf cmake
