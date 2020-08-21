# SLiSe: Framework to design rational filters for contour-based eigensolvers

We herein present a `Julia` implementation of `SLiSe`, a framework for the design of rational filters for contour-based eigensolvers.

The highlights include
- Unlimited _degrees of freedom_ within the design of rational filters,
- _Quick convergence_ with `SLiSe` filters through minimization of worst-case convergence rate,
- First toolkit to provide support for _box-constraints_, for fast convergence within contour-based eigensolvers with iterative linear solvers, and
- Fast, ready-to-use rational filters _included_.

## Installation
First, add this repository to `Julia` via
```jl
Pkg.clone("https://github.com/SimLabQuantumMaterials/SLiSeFilters.jl")
```

Afterwards, add the `L-BFGS-B` wrapper
```jl
Pkg.clone("https://github.com/jey/Lbfgsb.jl")
```
and run
```jl
Pkg.build("Lbfgsb")
```
to compile the original Fortran library.

## Usage
To obtain filters with small worst-case convergence rate, type into `Julia`:
```jl
import SLiSeFilters

# Choose poles per quadrant and gap size
p = 4
G = 0.95

# Starting filter
z, a   = SLiSeFilters.points( "zolo" , p );

# Starting weight function with shift
ws = SLiSeFilters.createWS(  [ sqrt(G),  sqrt(G)^-1, 1.4, 10 ],
                  BitVector( [ 0,    0,   1,  1 ] ) , # free vars within minimization
                             [ 1,  .01,  10, 20 ]   ,
                  BitVector( [ 0,    1,   1,  1 ] ) )   # free vars within minimization

# Compute filter with reduced worst-case convergence rate
z2, a2  = SLiSeFilters.wise( ws, z, a, G )
```
You find resulting, ready-to-use rational filters in the `docs/ ` folder, as well as Gauss-Legendre filters with reduced worst-case convergence rate.

A `Jupyter` notebook to demo the underlying `SLiSe` framework is available in the `docs/` folder.

## Related packages
A large set of benchmark eigenproblems for the [FEAST eigensolver](http://www.feast-solver.org/) can be obtained from the [SpectrumSlicingTestSuite.jl](https://github.com/SimLabQuantumMaterials/SpectrumSlicingTestSuite.jl) `Julia` package.

## References
- Jan Winkelmann and Edoardo Di Napoli. 2019. [Non-linear Least-Squares Optimization of Rational Filters for the Solution of Interior Hermitian Eigenvalue Problems.](https://doi.org/10.3389/fams.2019.00005) Frontiers in Applied Mathematics
  and Statistics, doi: 10.3389/fams.2019.00005.
- Eric Polizzi. 2009. [Density-matrix-based algorithm for solving eigenvalue problems.](https://doi.org/10.1103/PhysRevB.79.115112) Physical Review B, 79, 11, (Mar. 2009), 115112–115117. doi: 10.1103/PhysRevB.79.115112.
- Richard Byrd, Peihuang Lu, Jorge Nocedal, and Ciyou Zhu. 1995. [A limited memory algorithm for bound constrained optimization.](https://doi.org/10.1137/0916069) SIAM Journal on Scientific Computing, 16, 5, (Sept. 1995), 1190–1208. doi: 10.1137/0916069.
