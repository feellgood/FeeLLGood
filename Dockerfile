FROM ubuntu:22.04

LABEL description="This is an experimental image for running feeLLGood"
LABEL website="https://feellgood.neel.cnrs.fr/"

WORKDIR /root/src

# Copy only install-dependencies.sh
# This helps with layer caching
COPY ci-tests/install-dependencies.sh .

# Install dependencies
RUN apt-get update && apt-get install -y sudo python3 python3-pip
RUN bash install-dependencies.sh -u -d

# Copy source code
COPY . FeeLLGood

# Compile and install
RUN cd FeeLLGood && \
    sed -i 's/-march=native/-mavx2/' CMakeLists.txt && \
    cmake . && make -j $(nproc) && make install

# Remove build artefacts and source code
WORKDIR /
RUN rm -rf /root/src
