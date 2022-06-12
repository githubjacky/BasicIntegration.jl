include("../src/BasicIntegration.jl")
import .BasicIntegration: GHermite, GLaguerre, GLegendre, GQ
using QuadGK, Distributions, HCubature


function sepBeg(str; n=60)
    len = length(str)
    while n < len
        n += 5
    end
    u = floor(Int64, (n - len) / 2)
    a, b = floor(Int64, (3 / 4) * u), floor(Int64, (1 / 4) * u)
    k = n - (2(a + b) + len)
    while k > 1
        a += 1
        k -= 2
    end
    output = "#"^a * " "^b * str * " "^b * "#"^(a + k)
    println("#"^n)
    println(output)
    println("#"^n)
end


function sepEnd()
    println()
    println()
end
################################################################################
######################       one dimension problem      ########################
################################################################################
function errAnalysis1(g::Function, a, b; n=50)
    trueVal = quadgk(g, a, b)[1]
    # @show GHermite(g, a, b, n) - trueVal
    # @show GLaguerre(g, a, b, n) - trueVal
    @show GLegendre(g, a, b, n) - trueVal
    # @show GQ(g, b, n) - trueVal
end


# sepBeg("g(x)=exp(-x^2/3)*sqrt(1+x^2), A=[-Inf], b=[Inf]")
# g = x -> exp(-x^2 / 3) * sqrt(1 + x^2)
# a, b = -Inf, Inf
# errAnalysis1(g, a, b)
# sepEnd()
# #
# #
# sepBeg("g(x)=exp(-x^2/100)*(1+x^2), A=[0], b=[Inf]")
# g = x -> exp(-x^2 / 100) * (1 + x^2)
# a, b = 0, Inf
# errAnalysis1(g, a, b)
# sepEnd()
# #
# #
# sepBeg("g(x)=pdf(Normal(), x), A=[-Inf], b=[3]")
# g = x -> pdf(Normal(), x)
# a, b = -Inf, 3
# errAnalysis1(g, a, b)
# sepEnd()
# #
# #
sepBeg("g(x)=cdf(Normal(),2/3)^(-1)*pdf(Normal(2,3),x), A=[10], b=[10]")
g = x -> cdf(Normal(), 2 / 3)^(-1) * pdf(Normal(2, 3), x)
a, b = -10, 10
errAnalysis1(g, a, b)
sepEnd()


################################################################################
######################     multi-dimension problem      ########################
################################################################################
function errAnalysis2(g::Function, A::Vector{T}, B::Vector{P}; n=30) where {T,P<:Real}
    trueVal = hcubature(g, A, B, rtol=1e-8)[1]
    @show GHermite(g, A, B, n) - trueVal
    @show GLaguerre(g, A, B, n) - trueVal
    @show GLegendre(g, A, B, n) - trueVal
end

# sepBeg("multi dimension case(deviate), A=[0,0,0,0,0], B=[1,1,1,1,1]")
# g = x -> 1 / (x[1] + 1) + sqrt(x[2]) + 2 * (x[3]^2) + sqrt(2 * x[4]) + cbrt(x[5])
# A, B = [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]
# @btime errAnalysis2(g, A, B)
# sepEnd()


function errAnalysis3(g::Function, g1::Function, A::Vector{T}, B::Vector{P}; n=30) where {T,P<:Real}
    trueVal = quadgk(g1, A[1], B[1], rtol=1e-8)[1]
    # @show GHermite(g, A, B, n) - trueVal
    # @show GLaguerre(g, A, B, n) - trueVal
    @show GLegendre(g, A, B, n) - trueVal
end

# sepBeg("multi dimension case(deviate)")
# sepEnd()
#
# sepBeg("g(x)=exp(-x^2/3)*sqrt(1+x^2), A=[-Inf], b=[Inf]")
# g = x -> exp(-x[1]^2 / 3) * sqrt(1 + x[1]^2)
# g1 = x -> exp(-x^2 / 3) * sqrt(1 + x^2)
# a, b = -Inf, Inf
# errAnalysis3(g, g1, [a], [b])
# sepEnd()
# #
# #
# sepBeg("g(x)=exp(-x^2/100)*(1+x^2), A=[0], b=[Inf]")
# g = x -> exp(-x[1]^2 / 100) * (1 + x[1]^2)
# g1 = x -> exp(-x^2 / 100) * (1 + x^2)
# a, b = 0, Inf
# errAnalysis3(g, g1, [a], [b])
# sepEnd()
# #
# #
# sepBeg("g(x)=pdf(Normal(), x), A=[-Inf], b=[3]")
# g = x -> pdf(Normal(), x[1])
# g1 = x -> pdf(Normal(), x)
# a, b = -Inf, 3
# errAnalysis3(g, g1, [a], [b])
# sepEnd()
# #
# #
sepBeg("g(x)=cdf(Normal(),2/3)^(-1)*pdf(Normal(2,3),x), A=[10], b=[10]")
g = x -> cdf(Normal(), 2 / 3)^(-1) * pdf(Normal(2, 3), x[1])
g1 = x -> cdf(Normal(), 2 / 3)^(-1) * pdf(Normal(2, 3), x)
a, b = -10, 10
errAnalysis3(g, g1, [a], [b])
sepEnd()
