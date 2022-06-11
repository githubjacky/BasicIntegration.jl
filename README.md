# Basic Integraion
*pick the right method and don't need to worry about the domain*

### methodology
- Gaussian Quadrature(use [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) to generate the nodes and weights)
    - Gauss-Hermite
    - Gauss-Laguerre
    - Gauss-legendre
- Monte Carlo integration
    - vanilla Monte Carlo
    - Quassi Monte Carlo
    - Importance Sampling(can only use for one dimension problem just now)
        - using Halton sequences to do inverse transform sampling


### usage
```julia
# one dimension problem
g = x -> exp(-x^2 / 3) * sqrt(1 + x^2)
a, b = -Inf, Inf
GHermite(g, a, b, 30)  # the last parameter is the n-point Gauss Quadrature nodes and weights


# multidimension problem
g = x -> 1 / (x[1] + 1) + sqrt(x[2]) + 2 * (x[3]^2) + sqrt(2 * x[4]) + cbrt(x[5])
A, B = [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]
nodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863]
MCM(g, A, B, nodesNum, seed=1234)
quaMCM(g, A, B, nodesNum)
GLegendre(g, A, B, 30)
```

### implementation
```julia
function GHermite(g::Function, a::Real, b::Real, n::Int64) end
function GHermite(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real} end

function GLaguerre(g::Function, a::Real, b::Real, n::Int64) end
function GLaguerre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real} end

function GLegendre(g::Function, a::Real, b::Real, n::Int64) end
function GLegendre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real} end

#= provide another method to deal with the intergral: [-Inf, b],
    take advantage of both gausshermite and gausslaguerre =#
function GQ(g::Function, b::Real, n::Int64) end


function MCM(g::Function, a::Real, b::Real, nodesNum::Vector{Int64}; seed) end
function MCM(g::Function, A::Vector{T}, B::Vector{P}, nodesNum::Vector{Int64}; seed) where {T<:Real, P<:Real} end

function quaMCM(g::Function, a::Real, b::Real, nodesNum::Vector{Int64}) end
function quaMCM(g::Function, A::Vector{T}, B::Vector{P}, nodesNum::Vector{Int64}) where {T<:Real, P<:Real} end

function IShalton(g::Function, a::Real, b::Real, d::UnivariateDistribution, nodesNum::Vector{Int64}) end
```

### note
- newbie of julia programming but with enthusiasm
- still in progress...
