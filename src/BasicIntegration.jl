module BasicIntegration
export gausshermite, gausslaguerre, gausslegnedre, gaussmix
export montecarlo, quamontecarlo, importancesampling


import FastGaussQuadrature as FGQ
import HaltonSequences: HaltonPoint, Halton
import LinearAlgebra: dot
using Random, Distributions

include("GQutils.jl")
include("MCMutils.jl")


"""
Gauss Quadrature
"""
# one dimension case
function gausshermite(g, a::Real, b::Real, n::Integer=50)
    nodes, w = FGQ.gausshermite(n)
    invp(x) = exp(x^2)
    if a == -Inf && b == Inf
        f = x -> g(x) * invp(x)
    else
        x1, jcb1, x2, jcb2 = GHdt(a, b)
        f = t -> invp(t) * g(x1(x2(t))) * jcb1(x2(t)) * jcb2(t)
    end
    res = dot(f.(nodes), w)
    return res
end

gausshermite(g, d::Tuple{Any, Any}, n::Integer=30) = gausshermite(g, d[1], d[2], n)


# multidimension case
function gausshermite(g, A::Vector{T}, B::Vector{P}, n::Int) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    invp(x) = exp(x^2)
    nodes, w = FGQ.gausshermite(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausshermite nodes and weights
    f = MultiDim_GHdt(g, A, B, invp, dim)
    return dot(f.(Nodes), W)
end

gausshermite(g, d::Tuple{Vector{T}, Vector{T}}, n) where{T<:Real} = gausshermite(g, d[1], d[2], n)


function gausslaguerre(g::Function, a::Real, b::Real, n::Int64)  # one dimension case
    nodes, w = FGQ.gausslaguerre(n)
    invp(x) = exp(x)
    if a == 0 && b == Inf
        f = x -> invp(x) * g(x)
    else
        x1, jcb1, x2, jcb2 = GLadt(a, b)
        f = t -> invp(t) * g(x1(x2(t))) * jcb1(x2(t)) * jcb2(t)
    end
    return dot(f.(nodes), w)
end

gausslaguerre(g, d::Tuple{T, P}, n) where{T<:Real, P<:Real} = gausslaguerre(g, d[1], d[2], n)

# multidimension case
function gausslaguerre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    invp(x) = exp(x)
    nodes, w = FGQ.gausslaguerre(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausslaguerre nodes and weights
    f = MultiDim_GLadt(g, A, B, invp, dim)
    return dot(f.(Nodes), W)
end

gausslaguerre(g, d::Tuple{Vector{T}, Vector{T}}, n) where{T<:Real} = gausslaguerre(g, d[1], d[2], n)


function gausslegendre(g::Function, a::Real, b::Real, n::Int64)  # one dimension case
    nodes, w = FGQ.gausslegendre(n)
    if a == -1 && b == 1
        f = g
    else
        x1, jcb1 = GLdt(a, b)
        f = t -> g(x1(t)) * jcb1(t)
    end
    return dot(f.(nodes), w)
end

gausslegendre(g, d::Tuple{T, T}, n) where{T<:Real} = gausslegendre(g, d[1], d[2], n)

# multidimension case
function gausslegendre(g::Function, A::Vector{T}, B::Vector{P}, n::Int64) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    nodes, w = FGQ.gausslegendre(n)
    Nodes, W = nodesGen(n, dim, nodes, w)  # the algorithm to generate gausslegendre nodes and weights
    g = MultiDim_GLdt(g, A, B, dim)
    return dot(g.(Nodes), W)
end

gausslegendre(g, d::Tuple{Vector{T}, Vector{T}}, n) where{T<:Real} = gausslegendre(g, d[1], d[2], n)


#= provide another method to deal with the intergral: [-Inf, b],
    take advantage of both gausshermite and gausslaguerre =#
function gaussmix(g::Function, b::Real, n::Int64)
    nodes1, w1 = FGQ.gausshermite(n)
    invp1(x) = exp(x^2)
    f1 = x -> invp1(x) * g(x)
    nodes2, w2 = FGQ.gausslaguerre(n)
    invp2(x) = exp(x)
    f2 = b == 0 ? (x -> invp2(x) * g(x)) : (x -> invp2(x) * g(x+b))
    return dot(f1.(nodes1), w1) - dot(f2.(nodes2), w2)
end

gaussmix(g, d::Tuple{T, T}, n) where{T<:Real} = gaussmix(g, d[1], d[2], n)


function gaussmix(g::Function, A::Vector{T}, B::Vector{P}, n) where {T<:Real, P<:Real}
    length(A) == length(B) ? (dim = length(A)) : throw("dimension of domain inconsist")
    nodesArr = Vector{Vector{Float64}}(undef, 3) 
    wArr = Vector{Vector{Float64}}(undef, 3)
    nodesArr[1], wArr[1] = FGQ.gausshermite(n)
    nodesArr[2], wArr[2] = FGQ.gausslaguerre(n)
    nodesArr[3], wArr[3] = FGQ.gausslegendre(n)
    groupIndex = Vector{Int64}(undef, dim)
    g = addInvp(g, A, B, dim, groupIndex)
    Nodes, W = nodesGenGQ(n, dim, nodesArr, wArr, groupIndex)
    return dot(g.(Nodes), W)
end

gaussmix(g, d::Tuple{Vector{T}, Vector{T}}, n) where{T<:Real} = gaussmix(g, d[1], d[2], n)



"""
Monte Carlo Integration
"""
function montecarlo(g::Function, a::Real, b::Real; nodesNum::Vector{Int64}, seed)  # one dimension case
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
function montecarlo(g::Function, A::Vector{T}, B::Vector{P}, nodesNum::Vector{Int64}; seed) where {T<:Real, P<:Real}
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

# one dimension case
function quamontecarlo(g::Function, a::Real, b::Real; nodesNum::Vector{Int64})
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
function quamontecarlo(g::Function, A::Vector{T}, B::Vector{P}; nodesNum::Vector{Int64}) where {T<:Real, P<:Real}
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
function importancesampling(g::Function, a::Real, b::Real; dist::UnivariateDistribution, nodesNum::Vector{Int64})
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

end # end of BasicIntegration module