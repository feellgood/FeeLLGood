#!/bin/bash

# install-dependencies.sh: Install the dependencies required for
# building FeeLLGood.
#
# This is a support script for automated CI tests. It is not intended
# for end-users, and it contains a few calls to `sudo'. Use it at your
# own risk.
#
# Options:
#   -u  install the libraries required for unit tests
#   -d  install doxygen and graphviz

########################################################################
# The following few commands are not part of the documented installation
# procedure.

# Parse options.
while getopts "ud" OPT; do
    case "${OPT}" in
        u) unit_tests="true";;
        d) doxygen="true";;
        ?) exit 1;;
    esac
done

# Identify the OS.
eval $(grep -E '^(ID|PRETTY_NAME)=' /etc/os-release)
echo -e "\n*** Building on $PRETTY_NAME"

# Set shell options:
# -e: exit as soon as a command fails
# -v: print the lines being executed
set -ev

# Number of make jobs to run concurrently.
job_count=$(getconf _NPROCESSORS_ONLN)

########################################################################
# The following should match the documented installation procedure, save
# for the following options added to some commands:
# apt-get: -y, dnf: -y, make: -j, wget: -nv, unzip: -q

# Install required packages.
packages="unzip make cmake git"
if [ "$ID" = "rocky" ]; then
    packages="$packages tar wget gcc-c++ eigen3-devel tbb-devel yaml-cpp-devel duktape-devel"
    packages="$packages fftw-devel openblas-devel lapack-devel pkgconf"
    sudo dnf check-update -q || true
    sudo dnf install -y 'dnf-command(config-manager)'
    sudo dnf config-manager --set-enabled devel crb
    sudo dnf install -y epel-release
    if [ "$unit_tests" = "true" ]; then
        packages="$packages boost-devel"
    fi
    if [ "$doxygen" = "true" ]; then
        packages="$packages doxygen"
    fi
    sudo dnf install -y $packages
else  # Debian-like OS
    packages="$packages wget g++ libeigen3-dev libtbb-dev libyaml-cpp-dev duktape-dev"
    packages="$packages libann-dev libgmsh-dev"
    packages="$packages libfftw3-dev libopenblas-dev liblapack-dev pkg-config"
    sudo apt-get update -q
    if [ "$unit_tests" = "true" ]; then
        packages="$packages libboost-system-dev libboost-filesystem-dev libboost-test-dev"
    fi
    if [ "$doxygen" = "true" ]; then
        packages="$packages doxygen graphviz"
    fi
    sudo DEBIAN_FRONTEND=noninteractive apt-get install -y $packages
fi

# Download and build the libraries here.
mkdir -p ~/src
cd ~/src

# On Rocky Linux, download and install ANN and GMSH.
# On Debian and Ubuntu, they have already been installed with apt.
if [ "$ID" = "rocky" ]; then

    # Download and install ANN.
    rm -rf ann_1.1.2/
    if [ ! -f "ann_1.1.2.tar.gz" ]; then
        wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
    fi
    tar xzf ann_1.1.2.tar.gz
    cd ann_1.1.2/
    sed -i 's/CFLAGS =.* -O3/& -std=c++98/' Make-config
    make -C src -j $job_count linux-g++
    sudo cp lib/libANN.a /usr/local/lib/libann.a
    sudo cp --parents include/ANN/ANN.h /usr/local/
    cd ..

    # Download and install GMSH.
    if [ ! -f "gmsh-4.8.4-source.tgz" ]; then
        wget -nv https://gmsh.info/src/gmsh-4.8.4-source.tgz
    fi
    tar -xzf gmsh-4.8.4-source.tgz
    cd gmsh-4.8.4-source
    mkdir build
    cd build
    cmake -DDEFAULT=0 -DENABLE_BUILD_SHARED=1 ..
    make
    sudo make install/strip
    cd ../..
fi

# Download and install ScalFMM 3 (header only library). The --recursive option fetches its
# submodules (xtensor, xsimd, xtl, xtensor-blas, cpp_tools).
scalfmm_tag=V3.1.1
rm -rf ScalFMM/
git clone --depth 1 --branch $scalfmm_tag --recursive --shallow-submodules \
    https://gitlab.inria.fr/solverstack/ScalFMM.git
cmake -S ScalFMM -B ScalFMM/build -DCMAKE_BUILD_TYPE=Release \
    -Dscalfmm_BUILD_TOOLS=OFF -Dscalfmm_BUILD_EXAMPLES=OFF \
    -Dscalfmm_BUILD_CHECK=OFF -Dscalfmm_BUILD_UNITS=OFF
sudo cmake --install ScalFMM/build
