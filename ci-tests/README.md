# CI testing

This directory contains helper files for “continuous integration” tests,
i.e. for testing that feeLLGood builds and runs correctly on multiple
Linux distributions.

Automated tests are performed by [GitHub Actions][] on every push to the
repository using `ubuntu-latest` runners. These tests check that
feeLLGood builds and runs correctly, and that the unit tests pass. Tests
can also be run manually on a variety of Linux flavors. These manual
tests are meant to validate the published [installation procedure][]. We
check that feeLLGood builds and runs correctly, and that the doxygen
docs build.

Summary of tested items:

| test                     | automated | manual |
|--------------------------|-----------|--------|
| feeLLGood builds         |     ✔     |   ✔    |
| feeLLGood runs correctly |     ✔     |   ✔    |
| unit tests pass          |     ✔     |   ✘    |
| the documentation builds |     ✘     |   ✔    |

[GitHub Actions]: https://docs.github.com/en/actions
[installation procedure]: http://feellgood.neel.cnrs.fr/install.html

## Files

* `Vagrantfile`: used to spin up virtual machines for manual testing
* `install-dependencies.sh`: install the dependencies needed to build
  feeLLGood. This script accepts two options:
  * `-u`: install any extra dependencies needed to build the unit tests
  * `-d`: install doxygen, which is needed to build the docs
* `full_test.py`: minimal toy problem for checking that feeLLGood runs
  correctly

The script `full_test.py` uses the mesh `ellipsoid.msh` (163&nbsp;nodes
and 640&nbsp;elements), which is expected to be in the directory
`../examples`.

## Manual testing

Install [Vagrant][] and [VirtualBox][]:

```shell
curl -fsSL https://apt.releases.hashicorp.com/gpg | sudo apt-key add -
sudo apt-add-repository "deb [arch=amd64] https://apt.releases.hashicorp.com $(lsb_release -cs) main"
sudo apt install vagrant virtualbox
```

Then, in this directory, run

```shell
vagrant status
```

to see the list of supported virtual machine images. These are named
after the installed Linux distributions, e.g. `ubuntu_jammy` or
`rockylinux_9`. Next,

```shell
vagrant up machine
```

where ”machine” is one of the supported images. This will spin up a
virtual machine running the selected image, and provision it by running
`install-dependencies.sh -d` on it. Then

```shell
vagrant ssh machine
```

to open a shell on the VM.

If the image supports it, the parent of this directory is available
read-only within the VM as `host-FeeLLGood`. It can be copied in order
to test the working copy of the code, or cloned as a repository in order
to test a committed version. If `host-FeeLLGood` is not available, the
public repository should be cloned instead:

```shell
cp -a host-FeeLLGood FeeLLGood  # to test the working copy, or
git clone host-FeeLLGood FeeLLGood  # to test a commit, or
git clone https://github.com/feellgood/FeeLLGood.git  # if there is no host-FeeLLGood
```

The tests to run are:

```shell
cd FeeLLGood/
cmake .
make -j 2
sudo make install
doxygen
ci-tests/full_test.py
```

The VM can then be stopped by logging out of the ssh session and typing:

```shell
vagrant halt machine
```

[Vagrant]: https://www.vagrantup.com/
[VirtualBox]: https://www.virtualbox.org/
