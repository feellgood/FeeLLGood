# FEELLGOOD - micromagnetic solver.

FEELLGOOD is a micromagnetic solver using finite element technique to integrate Landau Lifshitz Gilbert equation. It computes the demagnetizing field using the so-called fast multipole algorithm.

It is developped by JC Toussaint & al.
This git repository is for devs, the code is being modified without any warranty it works. A public version can be found [here][]  

### Dependencies :
* C++ 11 and STL
* [ANN][] 1.1.2 version recommended
* [GMM][] 5.x version recommended (should be fine with 4.x too)
* [ScalFMM][] on branch release_1.5 tag v1.5.0 (version recommended)
* [GSL][] (optional)

[here]: http://feellgood.neel.cnrs.fr/
[ANN]: https://www.cs.umd.edu/~mount/ANN/
[GMM]: http://getfem.org/gmm.html/
[ScalFMM]: https://gitlab.inria.fr/solverstack/ScalFMM/
[GSL]: https://www.gnu.org/software/gsl/
