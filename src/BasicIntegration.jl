module BasicIntegration
export GHermite, GLaguerre, GLegendre, GQ, MCM, quaMCM, IShalton


import FastGaussQuadrature: gausshermite, gausslaguerre, gausslegendre
import HaltonSequences: HaltonPoint, Halton
import LinearAlgebra: dot
using Random, Distributions

include("GQutils.jl")
include("MCMutils.jl")


function GHermite(g::Function, a::Real, b::Real, n::Int64)  # one dimension case
    nodes, w = gausshermite(n)
    invp(x) = exp(x^2)
    if a == -Inf && b == Inf
        f = x -> g(x) * invp(x)
    else
        x1, jcb1, x2, jcb2 = GHdt(a, b)
        f = t -> invp(t) * g(x1(x2(t))) * jcb1(x2(t)) * jcb2(t)
    end
    return dot(f.(nodes), w)
end

# multidimension case
function GHermite(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    invp(x) = exp(x^2)
    nodes, w = gausshermite(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausshermite nodes and weights
    f = MultiDim_GHdt(g, A, B, invp, dim)
    return dot(f.(Nodes), W)
end


function GLaguerre(g::Function, a::Real, b::Real, n::Int64)  # one dimension case
    nodes, w = gausslaguerre(n)
    invp(x) = exp(x)
    if a == 0 && b == Inf
        f = x -> invp(x) * g(x)
    else
        x1, jcb1, x2, jcb2 = GLadt(a, b)
        f = t -> invp(t) * g(x1(x2(t))) * jcb1(x2(t)) * jcb2(t)
    end
    return dot(f.(nodes), w)
end

# multidimension case
function GLaguerre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    invp(x) = exp(x)
    nodes, w = gausslaguerre(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausshermite nodes and weights
    f = MultiDim_GLadt(g, A, B, invp, dim)
    return dot(f.(Nodes), W)
end


function GLegendre(g::Function, a::Real, b::Real, n::Int64)  # one dimension case
    nodes, w = gausslegendre(n)
    if a == -1 && b == 1
        f = g
    else
        x1, jcb1 = GLdt(a, b)
        f = t -> g(x1(t)) * jcb1(t)
    end
    return dot(f.(nodes), w)
end

# multidimension case
function GLegendre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    nodes, w = gausslegendre(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausshermite nodes and weights
    g = MultiDim_GLdt(g, A, B, dim)
    return dot(g.(Nodes), W)
end


#= provide another method to deal with the intergral: [-Inf, b],
    take advantage of both gausshermite and gausslaguerre =#
function GQ(g::Function, b::Real, n::Int64)
    nodes1, w1 = gausshermite(n)
    invp1(x) = exp(x^2)
    f1 = x -> invp1(x) * g(x)
    nodes2, w2 = gausslaguerre(n)
    invp2(x) = exp(x)
    f2 = b == 0 ? (x -> invp2(x) * g(x)) : (x -> invp2(x) * g(x+b))
    return dot(f1.(nodes1), w1) - dot(f2.(nodes2), w2)
end


function MCM(g::Function, a::Real, b::Real, nodesNum::Vector{Int64}; seed)  # one dimension case
    if a == 0 && b == 1
        f = g 
    else
        x1, jcb1 = MCMdt(a, b)
        f = t -> g(x1(t)) * jcb1(t)
    end
    n = length(nodesNum)
    res = Vector{Float64}(undef, n)
    for i = 1:n
        @inbounds res[i] = mean(f.(rand(Xoshiro(seed), nodesNum[i])))
    end
    return res
end

# multi dimension case
function MCM(g::Function, A::Vector{T}, B::Vector{P}, nodesNum::Vector{Int64}; seed) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : error("dimension of domain inconsist")
    n = length(nodesNum)
    f = multiDim_MCMdt(g, A, B, dim)
    res = Vector{Float64}(undef, n)
    for i = 1:n
        nodes = collect(eachrow(rand(Xoshiro(seed), nodesNum[i], dim)))
        @inbounds res[i] = mean(f.(nodes))
    end
    return res
end


function quaMCM(g::Function, a::Real, b::Real, nodesNum::Vector{Int64})  # one dimension case
    if a == 0 && b == 1
        f = g 
    else
        x1, jcb1 = MCMdt(a, b)
        f = t -> g(x1(t)) * jcb1(t)
    end
    n = length(nodesNum)
    res = Vector{Float64}(undef, n)
    for i = 1:n
        @inbounds res[i] = mean(f.(Halton(2, length=nodesNum[i])))
    end
    return res
end

# multi dimenstion case
function quaMCM(g::Function, A::Vector{T}, B::Vector{P}, nodesNum::Vector{Int64}) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : error("dimension of domain inconsist")
    n = length(nodesNum)
    f = multiDim_MCMdt(g, A, B, dim)
    res = Vector{Float64}(undef, n)
    for i = 1:n
        @inbounds res[i] = mean(f.(HaltonPoint(dim, length=nodesNum[i])))
    end
    return res
end


# one dimension case
function IShalton(g::Function, a::Real, b::Real, d::UnivariateDistribution, nodesNum::Vector{Int64})
    if a == 0 && b == 1
        f = g 
    else
        x1, jcb1 = MCMdt(a, b)
        f = t -> g(x1(t)) * jcb1(t)
    end
    F(x) = f(x) / pdf(d, x)
    n = length(nodesNum)
    res = Vector{Float64}(undef, n)
    for i = 1:n
        sample = quantile.(d, Halton(2, length=nodesNum[i]))
        @inbounds res[i] = mean(F.(sample))
    end
    return res
end

end # module
