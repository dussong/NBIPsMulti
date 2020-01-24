# NBIPsMulti

## Description

This Julia package extends the packages `cortner/NBodyIPs` and `cortner/IPFitting` to multicomponent systems. It generate Interatomic Potentials based on symmetry-adapted polynomials.

This is a *prototype* code.

## Installation notes

To install the package, first clone the repository:
```julia
] add https://github.com/dussong/NBIPsMulti.git
```
You may need to install additional packages.

To run the examples provided in the example folder, you need to activate the environment relative to the package. For this, go in the directory of the package and run
```julia
] activate .
```
Alternatively, you can run
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

For example, using bond-lengths PIPs with Morse coordinates `"exp( - 2.5 * (r/$r0-1))"` and a two-sided cosine cut-off `"(:cos2s, $(0.7*r0), $(0.88*r0), $(1.8*r0), $(rcut3))"`, we can define a descriptor
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
* `:AAAA, :AAAB, :AABC, :AABB, ABCD` (four-body)
We can then generate basis functions using `nbpolys`, e.g.,
```julia
#            descriptor, degree, atomic number of the species
B3 = nbpolys(BL3_AAA, 14, [6,6,6])
```
In practice, one would normally specify different cut-offs and space transforms
for different body-orders. Suppose these give descriptors `D2, D3aaa, D3aab`, with species types `:AA`, `:AAA` and `:AAB` then
a 4-body basis can be constructed providing the descriptor, polynomial degree and atomic number of the species involved (they have to match the descriptor species type)
```julia
B = [ nbpolys(D2, 14, [1,1]); nbpolys(D3aaa, 7, [1,1,1]); nbpolys(D3, 8, [1,1,6]) ]
```

For more details and more complex basis sets, see below.


### Step 3: Precompute a Lsq system

Once the dft dataset and the basis functions have been specified, the
least-squares system matrix can be assembled. This can be very time-consuming
for high body-orders, large basis sets and large data sets. Therefore this
matrix is stored in a block format that allows us to later re-use it in a variety
of different ways. This is done via
```julia
db = LsqDB(fname, configs, basis)
```
* The db is stored in two files: `fname_info.jld2` and `fname_kron.h5`. In
particular, `fname` is the path + name of the db, but without ending. E.g,
`"~/scratch/nbodyips/W5Bdeg12env"`.
* `configs` is a `Vector{Dat}`
* `basis` is a `Vector{<: AbstractCalculator}`
* The command `db = LsqDB(fname, configs, basis)` evaluates the basis functions,
e.g.,  `energy(basis[i], configs[j].at)` for all `i, j`, and stores these values
which make up the lsq matrix.

To reload a pre-computed lsq system, use `LsqDB(fname)`. To compute a lsq
system without storing it on disk, use `LsqDB("", configs, basis)`, i.e.,
pass an empty string as the filename.

### Step 4: Lsq fit, Analyse the fitting errors

The main function to call is
`lsqfit(db; kwargs...) -> IP, fitinfo`
The system is solved via (variants of) the QR factorisation. See `?lsqfit`
for details.

### Step 5: Usage

The output `IP` of `lsqfit` is a `JuLIP.AbstractCalculator` which supports
`energy, forces, virial, site_energies`. (todo: write more here, in
particular mention `fast`)


### More comments

- talk about the regularisation
- repulsive core.

there are two functions `filter_basis` and `filter_configs` that can be
used to choose a subset of the data and a subset of the basis. For example,
to take only 2B:
```
Ib2 = filter_basis(db, b -> (bodyorder(b) < 2))
```
See inline documentation for more details.


## Analysis

### Add fit information to a list of configurations

Suppose `configs::Vector{Dat}` is a list of configurations, then we can
add fitting error information by calling
```
add_fits!(myIP, configs, fitkey="myIP")
```
This will evaluate all observations stored in configs with the new IP and store
them in `configs[n].info[fitkey]["okey"]`. These observation values can then
be used to compute RMSE, produce scatter plots, etc.

This calculation can take a while. If `myIP` has just been fitted using `lsqfit`
then there is a quicker way to generate the fitting errors, but this is not
yet implemented. TODO: implement this!


## Hooks
