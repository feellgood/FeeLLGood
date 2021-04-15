#!/bin/bash

# install-dependencies.sh: Install the dependencies required for
# building FeeLLGood.
#
# This is a support script for the automated tests performed by GitHub
# Actions. It is not intended for end-users, and it contains a few calls
# to `sudo'. Use it at your own risk.

# Set shell options:
# -e: exit as soon as a command fails
# -v: print the lines being executed
set -ev

# Number of make jobs to run concurrently.
job_count=$(getconf _NPROCESSORS_ONLN)

# Install apt packages.
sudo apt-get update -q
sudo apt-get install -y libboost-dev \
    libboost-system-dev libboost-filesystem-dev libboost-test-dev

# Download and build the libraries here.
cd /tmp

# Install ANN.
wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/src/
make -j $job_count linux-g++
sudo cp ../lib/libANN.a /usr/local/lib/
sudo cp ../include/ANN/ANN.h /usr/local/include/
cd ../..

# Install exprtk.
wget -nv http://www.partow.net/downloads/exprtk.zip
unzip exprtk.zip
sudo cp exprtk/exprtk.hpp /usr/local/include/

# Patch, build and install ScalFMM.
wget -nv https://gitlab.inria.fr/solverstack/ScalFMM/-/archive/V1.5.1/ScalFMM-V1.5.1.tar.gz
tar xzf ScalFMM-V1.5.1.tar.gz
cd ScalFMM-V1.5.1/
sed -i 's/ROtation/Rotation/' Src/CMakeLists.txt
sed -i 's/~0x00LL/~0ULL/' Src/Containers/F{,Sub}Octree.hpp
sed -i 's/memcpy/if (nbParticles != 0) memcpy/' Src/Components/FBasicParticleContainer.hpp
cd Build
cmake ..
make -j $job_count
sudo make install
cd ../..

# Install GMM.
wget -nv http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.tar.gz
tar xzf gmm-5.4.tar.gz
cd gmm-5.4/
./configure
make -j $job_count
sudo make install
