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

# XXX: Test the behavior of `git show-ref' on Travis.
git show-ref
git show-ref HEAD
git show-ref --head HEAD
git show-ref --head HEAD --hash

# Stricter ref matching with `--verify'.
git show-ref --verify --head HEAD
git show-ref --verify --head HEAD --hash

# Force the build to fail early.
exit 1

# Install apt packages.
sudo apt-get update -q
sudo apt-get install -y cmake libboost-dev

# Download and build the libraries here. This is the grandparent of the
# current directory in a Travis build.
cd $HOME/build/

# Install ANN.
wget -nv https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xzf ann_1.1.2.tar.gz
cd ann_1.1.2/
make linux-g++
sudo cp lib/libANN.a /usr/local/lib/
sudo cp include/ANN/ANN.h /usr/local/include/
cd ..

# Install exprtk.
wget -nv http://www.partow.net/downloads/exprtk.zip
unzip exprtk.zip
sudo cp exprtk/exprtk.hpp /usr/local/include/

# Patch, build and install ScalFMM.
wget -nv https://gforge.inria.fr/frs/download.php/file/35369/SCALFMM-1.4-148.tar.gz
tar xzf SCALFMM-1.4-148.tar.gz
cd SCALFMM-1.4-148/
sed -i 's/ROtation/Rotation/' Src/CMakeLists.txt
mkdir Build
cd Build
cmake ..
make
sudo make install
cd ../..

# Install GMM.
wget -nv http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.3.tar.gz
tar xzf gmm-5.3.tar.gz
cd gmm-5.3/
./configure
make
sudo make install
