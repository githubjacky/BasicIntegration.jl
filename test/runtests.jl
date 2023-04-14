using Test, QuadGK, Distributions, HCubature
import BasicIntegration as BI

# using Revise, BenchmarkTools
# include("../src/BasicIntegration.jl")
# import .BasicIntegration as BI


################################################################################
###############################  utility function  #############################
################################################################################
function compare(g::Function, domain; n=50, gaussfunc=BI.gausshermite)
    res = gaussfunc(g, domain, n)
    trueval = quadgk(g, domain[1], domain[2])[1]
    @test  res ≈ trueval
end


################################################################################
######################       one dimension problem      ########################
################################################################################
@testset "Gauss Quadrature(univariate)" begin
    compare(
        x -> exp(-x^2 / 3) * sqrt(1 + x^2),
        (-Inf, Inf),
        gaussfunc=BI.gausshermite
    )
    # compare(
    #     x -> exp(-x^2 / 100) * (1 + x^2),
    #     (0, Inf),
    #     gaussfunc=BI.gausslaguerre
    # )
    # compare(
    #     x -> pdf(Normal(), x),
    #     (-Inf, 3),
    #     gaussfunc=BI.gaussmix
    # )
    # compare(
    #     x -> cdf(Normal(), 2 / 3)^(-1) * pdf(Normal(2, 3), x),
    #     (-10, 10),
    #     gaussfunc=BI.gausslegendre
    # )
end