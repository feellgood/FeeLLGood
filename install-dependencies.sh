#!/bin/bash

# install-dependencies.sh: Install the dependencies required for
# building FeeLLGood.
#
# This is a support script for the automated tests performed by GitHub
# Actions. It is not intended for end-users, and it contains a few calls
# to `sudo'. Use it at your own risk.

########################################################################
# The following few commands are not part of the documented installation
# procedure.

# Set shell options:
# -e: exit as soon as a command fails
# -v: print the lines being executed
set -ev

# Number of make jobs to run concurrently.
job_count=$(getconf _NPROCESSORS_ONLN)

########################################################################
# The following should match the documented installation procedure, save
# for the following options added to some commands:
# apt-get: -y, make: -j, wget: -nv, unzip: -q

# Install required packages.
sudo apt-get update -q
sudo apt-get install -y libboost-dev \
    libboost-system-dev libboost-filesystem-dev libboost-test-dev

# Download and build the libraries here.
cd /tmp

# Download and install ANN.
wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/
make -j $job_count linux-g++
sudo cp lib/libANN.a /usr/local/lib/
sudo cp include/ANN/ANN.h /usr/local/include/
cd ..

# Download and install exprtk.
# Use the version tagged "g20210303" by FreeBSD ports.
# C.f. https://cgit.freebsd.org/ports/tree/math/exprtk/Makefile
exprtk_sha1=ca5c577917646ddba3f71ce6d5dd7d01f351ee80
wget -nv https://github.com/ArashPartow/exprtk/archive/$exprtk_sha1.zip
mv $exprtk_sha1.zip exprtk-$exprtk_sha1.zip
unzip -q exprtk-$exprtk_sha1.zip
sudo cp exprtk-$exprtk_sha1/exprtk.hpp /usr/local/include/

# Download, patch and install ScalFMM.
# Use the tip of the branch maintenance/scalfmm-1.5 as of 2021-09-20.
scalfmm_sha1=22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc
wget -nv https://gitlab.inria.fr/solverstack/ScalFMM/-/archive/$scalfmm_sha1/ScalFMM-$scalfmm_sha1.tar.gz
tar xzf ScalFMM-$scalfmm_sha1.tar.gz
cd ScalFMM-$scalfmm_sha1/
sed -i 's/memcpy/if (nbParticles != 0) memcpy/' Src/Components/FBasicParticleContainer.hpp
cd Build
cmake ..
make -j $job_count
sudo make install
cd ../..

# Download and install GMM.
wget -nv http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.tar.gz
tar xzf gmm-5.4.tar.gz
cd gmm-5.4/
./configure
make -j $job_count
sudo make install
