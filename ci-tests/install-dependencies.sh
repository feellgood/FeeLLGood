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
    packages="$packages wget gcc-c++ eigen3-devel tbb-devel yaml-cpp-devel"
    sudo dnf check-update -q || true
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
    if [ "$ID" = "debian" ]; then  # needed for libmkl-dev on Debian
        sudo apt-get install -y software-properties-common
        sudo apt-add-repository -y non-free
    fi
    packages="$packages wget g++ libeigen3-dev libmkl-dev libtbb-dev libyaml-cpp-dev duktape-dev"
    packages="$packages libgmsh-dev"
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

# Download and install ANN.
rm -rf ann_1.1.2/
if [ ! -f "ann_1.1.2.tar.gz" ]; then
    wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
fi
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/
sed -i 's/CFLAGS =.* -O3/& -std=c++98/' Make-config
make -j $job_count linux-g++
sudo cp lib/libANN.a /usr/local/lib/
sudo cp include/ANN/ANN.h /usr/local/include/
cd ..

# On Rocky Linux, download and install Duktape and GMSH.
# On Debian and Ubuntu, it has already been installed with apt.
if [ "$ID" = "rocky" ]; then
    if [ ! -f "duktape-2.7.0.tar.xz" ]; then
        wget -nv https://duktape.org/duktape-2.7.0.tar.xz
    fi
    tar -xJf duktape-2.7.0.tar.xz
    cd duktape-2.7.0/src
    gcc -O2 -c duktape.c
    ar rcs libduktape.a duktape.o
    sudo cp libduktape.a /usr/local/lib/
    sudo cp duktape.h duk_config.h /usr/local/include/
    cd ../..

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

# Download, patch and install ScalFMM.
# Use the tip of the branch maintenance/scalfmm-1.5 as of 2021-09-20.
scalfmm_sha1=22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc
rm -rf ScalFMM-$scalfmm_sha1/
if [ ! -f "ScalFMM-$scalfmm_sha1.tar.gz" ]; then
    wget -nv https://gitlab.inria.fr/solverstack/ScalFMM/-/archive/$scalfmm_sha1/ScalFMM-$scalfmm_sha1.tar.gz
fi
tar xzf ScalFMM-$scalfmm_sha1.tar.gz
cd ScalFMM-$scalfmm_sha1/
sed -i 's/memcpy/if (nbParticles != 0) memcpy/' Src/Components/FBasicParticleContainer.hpp
sed -i 's/OPENMP_CXX_FOUND/OPENMP_FOUND OR OPENMP_CXX_FOUND/' CMakeLists.txt
cd Build
cmake ..
make -j $job_count
sudo make install
