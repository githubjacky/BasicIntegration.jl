include("../src/BasicIntegration.jl")
import .BasicIntegration: MCM, quaMCM, IShalton
using QuadGK, Distributions


################################################################################
######################       one dimension problem      ########################
################################################################################
# nodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863]
# g = x -> exp(-x^2 / 3) * sqrt(1 + x^2)
# a, b = -Inf, Inf
# @show MCM(g, a, b, nodesNum, seed=1234)
# @show quaMCM(g, a, b, nodesNum)
# @show quadgk(g, a, b)[1]


g(x, c=1e-9, k=2) = c * x^(-k-1) * (1-x)^(k+1)
a, b = 1e-5, 1
nodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863]
d = truncated(Normal(1e-3, 1e-3), 0, 1)
@show IShalton(g, a, b, d, nodesNum)
@show quadgk(g, a, b)[1]

################################################################################
######################     multi-dimension problem      ########################
################################################################################
# nodesNum = [30^5]
# g = x -> 1 / (x[1] + 1) + sqrt(x[2]) + 2 * (x[3]^2) + sqrt(2 * x[4]) + cbrt(x[5])
# A, B = [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]
# @show MCM(g, A, B, nodesNum, seed=1234)
# @show quaMCM(g, A, B, nodesNum)
# @show hcubature(g, A, B, rtol=1e-8)[1]


# g = x -> sqrt(x[1]^2+x[2]^2+x[3]^2) * exp(-(x[1]^2+x[2]^2+x[3]^2))
# A, B = [-Inf, -Inf, -Inf], [Inf, Inf, Inf]
# nodesNum = [2^19-1, 2^25-1]
# res = MCM(g, A, B, nodesNum, seed=1234)

