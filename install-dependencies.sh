#!/bin/bash

# install-dependencies.sh: Install the dependencies required for
# building FeeLLGood.
#
# This is a support script for the automated build tests performed by
# Travis CI (travis-ci.org). It is not intended for end-users, and it
# contains a few calls to `sudo'. Use it at your own risk.

# Set shell options:
# -e: exit as soon as a command fails
# -v: print the lines being executed
set -ev

# Install apt packages.
sudo apt-get update -q
sudo apt-get install -y cmake libboost-dev

# Install the libraries in sibling directories of FeeLLGood.
cd ..

# Install ANN.
wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/
make linux-g++
sudo cp lib/libANN.a /usr/local/lib/
cd ..

# Install exprtk.
wget -nv http://www.partow.net/downloads/exprtk.zip
unzip exprtk.zip

# Install ScalFMM.
wget -nv https://gforge.inria.fr/frs/download.php/file/35369/SCALFMM-1.4-148.tar.gz
tar xzf SCALFMM-1.4-148.tar.gz
mkdir SCALFMM-1.4-148/Build
cd SCALFMM-1.4-148/Build
cmake ..
make
sudo make install
cd ../..

# Install GMM.
wget -nv http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.3.tar.gz
tar xzf gmm-5.3.tar.gz
