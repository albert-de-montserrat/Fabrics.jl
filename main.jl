import Pkg; Pkg.activate(".")
using Fabrics

# Example: Weak Inclusion
Δη, ϕ = 1e-3, 0.1
α_bulk = 30.5
p = parameterization(Δη, ϕ)
@show r1, r2 = fabric_shape(10, 1, p)
@show α = fabric_orientation(α_bulk, p)
@show η_cutoff = weakening(Δη, ϕ)

# Example: Strong Inclusion
Δη, ϕ = 1e1, 0.2
α_bulk = 30.5
p = parameterization(Δη, ϕ)
@show r1, r2 = fabric_shape(10, 1, p)
@show α = fabric_orientation(α_bulk, p)