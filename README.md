# FeeLLGood – A micromagnetic solver ![Build Status](https://github.com/feellgood/FeeLLGood/actions/workflows/tests.yml/badge.svg)

FEELLGOOD is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation. It computes the demagnetizing field using the so-called fast multipole algorithm.

It is developped by JC Toussaint & al.
The code is being modified without any warranty it works. A dedicated website can be found [here][]  

### Dependencies

* C++17 and the STL
* [TBB][]
* [yaml-cpp][]
* [ANN][] 1.1.2
* [Duktape][] 2.7.0
* [ScalFMM][] revision [22b9e4f6cf][ScalFMM-rev] (it should also work with V1.5.1)
* [Eigen][] ≥ 3.3

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
[ScalFMM-rev]: https://gitlab.inria.fr/solverstack/ScalFMM/-/archive/22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc/ScalFMM-22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc.tar.gz
[Eigen]: https://eigen.tuxfamily.org/
