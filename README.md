# FeeLLGood – A micromagnetic solver ![Build Status](https://github.com/feellgood/FeeLLGood/actions/workflows/tests.yml/badge.svg)

FEELLGOOD is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation. It computes the demagnetizing field using the so-called fast multipole algorithm.

It is developped by JC Toussaint & al.
The code is being modified without any warranty it works. A dedicated website can be found [here][]  
We recommand to use branch master, other branches are experimental or work in progress.

### Dependencies

* C++17 and the STL
* [TBB][]
* [yaml-cpp][]
* [ANN][] 1.1.2
* [Duktape][] 2.7.0
* [ScalFMM][] 3 (tested with V3.1.1), which requires [FFTW][], BLAS and LAPACK
* [Eigen][] ≥ 3.3
* [GMSH][] ≥ 4.8

ScalFMM 3 is a header only library, its submodules must be fetched with it:

```shell
git clone --branch V3.1.1 --recursive https://gitlab.inria.fr/solverstack/ScalFMM.git
cmake -S ScalFMM -B ScalFMM/build -DCMAKE_BUILD_TYPE=Release -Dscalfmm_BUILD_TOOLS=OFF \
      -Dscalfmm_BUILD_EXAMPLES=OFF -Dscalfmm_BUILD_CHECK=OFF -Dscalfmm_BUILD_UNITS=OFF
sudo cmake --install ScalFMM/build
```

If ScalFMM is installed elsewhere (`-DCMAKE_INSTALL_PREFIX=<prefix>`), configure feeLLGood with
`cmake . -Dscalfmm_DIR=<prefix>/lib/cmake/scalfmm`.

Intel MKL is not used by default: when Eigen uses MKL, the MKL cblas header conflicts with the one
of xflens, used by ScalFMM 3. feeLLGood selects OpenBLAS, else the generic BLAS and LAPACK; another
vendor can be chosen with `-DBLA_VENDOR=...`.

The parameters of the fast multipole method (interpolation order, tree height, group size) can be
set in the section `demagnetizing_field_solver` of the settings, see `feellgood --print-defaults`.

### License

Copyright (C) 2012-2023  Jean-Christophe Toussaint, with contributions by F. Alouges, D. Gusakova, S. Jamet, M. Struma, C. Thirion and E. Bonet.

FeeLLGood is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Additional permission under GNU GPL version 3 section 7: If you modify this Program, or any covered work, by linking or combining it with the Intel® MKL library (or a modified version of that library), containing parts covered by the terms of Intel Simplified Software License, the licensors of this Program grant you additional permission to convey the resulting work.

The libraries used by feeLLGood are distributed under different licenses, and this is documented in their respective Web sites.

[here]: https://feellgood.neel.cnrs.fr/
[TBB]: https://www.threadingbuildingblocks.org/
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[ANN]: https://www.cs.umd.edu/~mount/ANN/
[Duktape]: https://duktape.org/
[ScalFMM]: https://gitlab.inria.fr/solverstack/ScalFMM/
[FFTW]: https://www.fftw.org/
[Eigen]: https://eigen.tuxfamily.org/
[GMSH]: http://gmsh.info/
