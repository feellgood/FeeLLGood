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

[here]: https://feellgood.neel.cnrs.fr/
[TBB]: https://www.threadingbuildingblocks.org/
[yaml-cpp]: https://github.com/jbeder/yaml-cpp
[ANN]: https://www.cs.umd.edu/~mount/ANN/
[Duktape]: https://duktape.org/
[ScalFMM]: https://gitlab.inria.fr/solverstack/ScalFMM/
[ScalFMM-rev]: https://gitlab.inria.fr/solverstack/ScalFMM/-/archive/22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc/ScalFMM-22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc.tar.gz
[Eigen]: https://eigen.tuxfamily.org/
