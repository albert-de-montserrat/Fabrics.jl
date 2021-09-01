# Fabrics.jl

## Installation

```julia
using Pkg
Pkg.add("https://github.com/albert-de-montserrat/Fabrics.jl.git")
```

## Example

First we define the viscosity contrast and volume fraction of the inclusion phase
```julia
Δη, ϕ = 1e-3, 0.1
```
Obtain parameterization object
```julia
p = parameterization(Δη, ϕ)
```
We need the longest and intermediate semi-axes of the bulk Finite Strain Ellipsoid (FSE)
```julia
a1, a2 = 10, 1
```
Now we can calculate the fabric shape and viscosity cut-off
```julia
r1, r2 = fabric_shape(10, 1, p)
η_cutoff = weakening(Δη, ϕ)
```
We can also calculate the orientation of the fabric, given the orientation of the first axis of the bulk FSE
```julia
α_bulk = 30.5
α = fabric_orientation(α_bulk, p)
```
