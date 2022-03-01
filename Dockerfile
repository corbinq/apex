FROM ubuntu:20.04 as base

LABEL org.label-schema.name="apex"
LABEL org.label-schema.description="Toolkit for analysis of molecular quantitative trait loci (xQTLs)"
LABEL org.label-schema.url="https://github.com/corbinq/apex"
LABEL org.label-schema.usage="https://corbinq.github.io/apex/"
LABEL org.label-schema.vcs-url="https://github.com/corbinq/apex"
LABEL org.label-schema.schema-version="1.0"

# Install required packages for apex to install.
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  cmake \
  python3-pip \
  zlib1g-dev \
  liblzma-dev \
  libbz2-dev \
  locales \
  git \
  pkg-config \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Install python dependencies
RUN pip3 install --upgrade pip && \
    pip3 install cget

# Create a group and user to execute as, then drop root
ARG UID
ARG GID
RUN \
  if [ -n "$GID" ]; then \
    addgroup --gid $GID apex; \
  else \
    addgroup apex; \
  fi && \
  if [ -n "$UID" ]; then \
    adduser --gecos "User for running apex as non-root" --shell /bin/bash --disabled-password --uid $UID --ingroup apex apex; \
  else \
    adduser --gecos "User for running apex as non-root" --shell /bin/bash --disabled-password --ingroup apex apex; \
  fi

WORKDIR /home/apex
USER apex

# Install cpp dependencies
COPY --chown=apex:apex requirements.txt /home/apex/requirements.txt
COPY --chown=apex:apex *.cmake /home/apex/
ARG CMAKE_BUILD_PARALLEL_LEVEL
ARG MAKEFLAGS
RUN cget install -f requirements.txt

# Next stage: compile
FROM base as compile

# Copy source
COPY --chown=apex:apex CMakeLists.txt /home/apex/CMakeLists.txt
COPY --chown=apex:apex ./src /home/apex/src
COPY --chown=apex:apex ./tests /home/apex/tests

# Run compile
ENV CGET_PREFIX="/home/apex/cget"
ENV INSTALL_PREFIX="/home/apex/cget"
RUN \
  mkdir build \
  && cd build \
  && cmake .. \
    -DCMAKE_TOOLCHAIN_FILE=${CGET_PREFIX}/cget/cget.cmake \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
    -DCMAKE_BUILD_TYPE=Release
RUN cd build && cmake --build . --target apex
RUN cd build && cmake --build . --target apex-bin
RUN cd build && cmake --build . --target tests

# Run test cases
FROM compile as test
RUN cd tests && ../build/tests/tests

# Frequently changing metadata here to avoid cache misses
ARG BUILD_DATE
ARG GIT_SHA
ARG APEX_VERSION

LABEL org.label-schema.version=$apex_VERSION \
      org.label-schema.vcs-ref=$GIT_SHA \
      org.label-schema.build-date=$BUILD_DATE

# Set the default stage to be the base files + compiled binaries + test cases.
FROM test
