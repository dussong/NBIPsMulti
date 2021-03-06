# NBIPsMulti

## Description

This Julia package extends the packages `cortner/NBodyIPs.jl` and `cortner/IPFitting.jl` to multicomponent systems. It generate Interatomic Potentials based on symmetry-adapted polynomials.

This is an *experimental* code.

## Installation notes

To install the package, first clone the repository:
```julia
] dev https://github.com/dussong/NBIPsMulti.git
```
You may need to install additional packages. In particular, the package manager `https://github.com/JuliaMolSim/MolSim` may be useful.

To run the examples provided in the example folder, one first needs to activate the environment relative to the package. For this, run in the directory of the package
```julia
] activate .
```
Alternatively, one can run
```julia
] activate dir/NBodyIPsMulti
```
where `dir` is the directory where the package is cloned.


## Basic Usage

### Step 1: Import data/observations

Data importation is done with the `IPFitting` package.

### Step 2: Generate a basis

A basis is defined by
* choice of bond-length or bond-angle PIPs.
* space transform
* choice of cut-off
* body-order
* species type (main difference from the single species case)
* polynomial degree

For example, using bond-lengths PIPs with Morse coordinates `"exp( - 2.5 * (r/$r0-1))"` and a two-sided cosine cut-off `"(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))"`, one can define a descriptor
```julia
r0 = 3*round(rnn(:C),digits=2)
rcut2 = 2.8 * r0
rcut3 = 2.3 * r0
rcut4 = 1.9 * r0
BL3_AAA = MultiDesc("exp( - 2.5 * (r/$r0-1))", "(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))",Val(:AAA))
```
The species type `:AAA` means that the three particles are the same. The possible species types are
* `:AA` (two-body)
* `:AAA, :AAB, :ABC` (three-body)
* `:AAAA, :AAAB, :AABC, :AABB, :ABCD` (four-body)
One can then generate basis functions using `nbpolys`, providing the descriptor, polynomial degree and atomic number of the species involved (they have to match the descriptor species type), e.g.,
```julia
#            descriptor, degree, atomic number of the species
B3 = nbpolys(BL3_AAA, 14, [6,6,6])
```
In practice, one would normally specify different cut-offs and space transforms
for different body-orders. Suppose these give descriptors `D2, D3aaa, D3aab`, with species types `:AA`, `:AAA` and `:AAB` then
a 4-body basis can be constructed with
```julia
B = [ nbpolys(D2, 14, [1,1]); nbpolys(D3aaa, 7, [1,1,1]); nbpolys(D3, 8, [1,1,6]) ]
```

### Step 3: Precompute a Lsq system

As in `IPFitting.jl` the assembly of the least-square system is done via
```julia
db = LsqDB(fname, configs, basis)
```
see `IPFitting.jl` for more documentation.

### Step 4: Lsq fit, Analyse the fitting errors

The main function to call is
`lsqfit(db; kwargs...) -> IP, fitinfo`
The system is solved via (variants of) the QR factorisation. See `?lsqfit`
for details.

### Step 5: Usage

The output `IP` of `lsqfit` is a `JuLIP.AbstractCalculator` which supports
`energy, forces, site_energies`.

### More comments

There are some more capabilities, namely:
* Regularization for the resolution of the least-squares system. This is illustrated in the example `fit_with_regularisation.jl`
* Addition of a repulsive core for small interatomic distances. This is presented in the example `fit_repulsive_core.jl`
* Environment dependence example in `fit_with_environment.jl` (experimental).
