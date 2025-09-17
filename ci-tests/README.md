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
[installation procedure]: https://feellgood.neel.cnrs.fr/install.html

## Files

* `test-vm`: shell script for managing the virtual machines used for
  manual testing
* `install-dependencies.sh`: install the dependencies needed to build
  feeLLGood. This script accepts two options:
  * `-u`: install any extra dependencies needed to build the unit tests
  * `-d`: install doxygen, which is needed to build the docs
* `full_test.py`: minimal toy problem for checking that feeLLGood runs
  correctly

The script `full_test.py` uses the mesh `ellipsoid.msh` (167&nbsp;nodes
and 773&nbsp;elements), which is expected to be in the directory
`../examples`.

## Manual testing

Manual tests are done on VMs (virtual machines). In order to run these
VMs, install [Incus][], [QEMU][] and genisoimage:

```shell
sudo apt install incus qemu-system-x86 genisoimage
```

Then, in this directory, run

```shell
./test-vm list
```

to see the list of supported VMs. These are named after the installed
Linux distributions, e.g. “noble” for Ubuntu 24.04 LTS (Noble Numbat).
Each machine has a state that can be either “NOT CREATED”, “STOPPED” or
“RUNNING”. A machine that is not “RUNNING” can be started with:

```shell
machine=...
./test-vm start $machine
```

with `...` replaced by the machine name. If the machine was “NOT
CREATED”, this will download a machine image and create the VM from it,
which may take a few minutes.

Now, the feeLLGood dependencies can be installed by running

```shell
./test-vm provision $machine
```

If a directory named “dl” exists as a sibling of the test-vm script, it
will be used as a download cache: its content will be copied to
`/home/user/src` on the VM. Putting there the tarballs of ANN and
ScalFMM will prevent them from being downloaded by
install-dependencies.sh running inside the VM.

Finally,

```shell
./test-vm shell $machine
```

will open a shell on the VM.

On this machine, you can clone the feeLLGood repository:

```shell
cd src
git clone https://github.com/feellgood/FeeLLGood.git
```

The tests to run are:

```shell
cd FeeLLGood/
cmake .
make -j 2
ci-tests/full_test.py
sudo make install
doxygen > /dev/null
```

The VM can then be stopped by logging out of the shell session and typing:

```shell
./test-vm stop $machine
```

A short usage summary of the `test-vm` script can be displayed with:

```shell
./test-vm help
```

[Incus]: https://linuxcontainers.org/incus/
[QEMU]: https://www.qemu.org/
