module Fabrics

# exports 
export ParameterizationWeak
export ParameterizationStrong
export Parameterization
export parameterization
export fabric_shape
export fabric_orientation
export weakening
export Weakening

# Types and structures
abstract type Parameterization end

struct ParameterizationWeak{T} <: Parameterization
    ζ::T # Fabric shape
    α::T # Fabric orientation
end

struct Weakening{T} <: Parameterization
    ζ::T
    ξ::T
    χ::T
    θ::T
    ψ::T
    λ::T
end

struct ParameterizationStrong{T} <: Parameterization
    ζ::T # Fabric shape
    ξ::T # Fabric shape
    χ::T # Fabric shape
    θ::T # Fabric shape
    α::T # Fabric orientation
end

# Fabric shape (weak inclusions)
function fabric_shape(a1, a2, p::ParameterizationWeak{T}) where T
    A = log10(a1/a2) # Bulk semi axes ratios
    r₁ = p.ζ*A  # inclusion  log10(a1/a2)
    r₂ = A # inclusion  log10(a2/a3)
    return r₁, r₂ 
end

# Fabric shape (weak inclusions)
function fabric_shape(a1, a2, p::NamedTuple) 
    A = log10(a1/a2) # Bulk semi axes ratios
    R1, R2 = p.R1, p.R2
    r₁ = R1.ζ + R1.ξ*A + R1.χ*A^2 + R1.θ*A^3 # inclusion  log10(a1/a2)
    r₂ = R2.ζ + R2.ξ*A + R2.χ*A^2 + R2.θ*A^3 # inclusion  log10(a2/a3)
    return r₁, r₂ 
end

# Fabric orientation: weak inclusions
fabric_orientation(α,  p::T) where T<:Parameterization = p.α*α

# Fabric orientation: strong inclusions
fabric_orientation(α,  p::NamedTuple) = p.R1.α*α

# Proportionality constanst of shape and orientation for 2-phase Newtonian aggregate 
# with weak inclusion. Numbers from de Montserrat et al. 2021, JGR, Table 1 
parameterization(Δη::T, ϕ::T) where T<:Real = _barrier(Δη, ϕ)

_barrier(Δη, ϕ) = parameterization(Val{Δη}(), Val{ϕ}())

parameterization(::Val{1e-3}, ::Val{0.1}) = ParameterizationWeak(
    0.687903 ,
    1.000936
)

parameterization(::Val{1e-3}, ::Val{0.2}) = ParameterizationWeak(
    0.732392 ,
    0.933471
)

parameterization(::Val{1e-3}, ::Val{0.3}) = ParameterizationWeak(
    0.709184 ,
    1.095165
)

parameterization(::Val{1e-2}, ::Val{0.1}) = ParameterizationWeak(
    0.715399 ,
    1.
)

parameterization(::Val{1e-2}, ::Val{0.2}) = ParameterizationWeak(
    0.774659 ,
    1.002619
)

parameterization(::Val{1e-2}, ::Val{0.3}) = ParameterizationWeak(
    0.704054 ,
    0.945960
)

parameterization(::Val{1e-1}, ::Val{0.1}) = ParameterizationWeak(
    0.793576 ,
    1.07913
)

parameterization(::Val{1e-1}, ::Val{0.3}) = ParameterizationWeak(
    0.752338 ,
    0.983557
)

# Proportionality constanst of shape and orientation for 2-phase Newtonian aggregate 
# with strong inclusions. Numbers from de Montserrat et al. 2021, JGR, Table 2 
function parameterization(::Val{1e1}, ::Val{0.1}) 
    R1 = ParameterizationStrong(
        -0.006131,  0.458978, -0.125306, -0.083918, 1.0
    )
    R2 = ParameterizationStrong(
        0.000512, -0.319324 , 0.185224 , -0.306101, 1.0
    )
    return (R1=R1, R2=R2)
end

function parameterization(::Val{1e1}, ::Val{0.2}) 
    R1 = ParameterizationStrong(
        -0.000913 , 0.389964,-0.056822,0.1010262, 1.020180
    )
    R2 = ParameterizationStrong(
        0.002493 , 0.296955, 0.271349, -0.2545942, 1.020180
    )
    return (R1=R1, R2=R2)
end

function parameterization(::Val{1e1}, ::Val{0.3}) 
    R1 = ParameterizationStrong(
        0.000295 , 0.372009,-0.070759, 0.2408665, 1.043904
    )
    R2 = ParameterizationStrong(
        0.005547 , 0.266632, 0.271325, -0.1762204, 1.043904
    )
    return (R1=R1, R2=R2)
end

# Viscosity cut-off for 2-phase Newtonian aggregates with weak inclusion. 
# Numbers from de Montserrat et al. 2021, JGR, Table 1 
weakening(Δη::T, ϕ::T) where T<:Real = weakening(Δη, _barrier_w(ϕ))

_barrier_w(ϕ) = weakening(Val{ϕ}())

weakening(::Val{0.1}) = Weakening(
    0.903221, 1.20365, -1.148182, 0.3735, 0.698558, 0.6619,
)

weakening(::Val{0.2}) = Weakening(
    43.646559, -42.907958, 0.0, 0.261095, 0.984129, 1.0,
)

weakening(::Val{0.3}) = Weakening(
    40.628095, -39.823767, 0.0, 0.196188, 0.985677, 1.0,
)

weakening(Δη, w::Weakening) =  w.ζ*Δη^w.ψ + w.ξ*Δη^w.λ + w.χ*Δη + w.θ

end # module