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
        - use [Halton sequences](https://github.com/tobydriscoll/HaltonSequences.jl) to implment
        - since the prime I choose is 2 so it's better to pick nodes$=2^n-1$
    - Importance Sampling
        - input should include distribution
        - using Halton sequences to do inverse transform sampling


### function
- `GHermite()`, `GLaguerre()`, `GLegendre()`, `GQ()`
    - `GQ()` takes advantage of gausshermite and gausslaguerre to deal with the problem with domain [-Inf, 0] 
        - only support univariable problem
    - all the others can estimate multivariale problem with the vector-like input
- `MCM()`, `quaMCM()`, `IShalton()`
    - `IShalton()` can only be used for univariable problem just now


### usage
```julia
# univariable
g = x -> exp(-x^2 / 3) * sqrt(1 + x^2)
a, b = -Inf, Inf  # domain
GHermite(g, a, b, 30)  # the last parameter is the n-point Gauss Quadrature nodes and weights
GLaguerre(g, a, b, 30) 

m(x, c=1e-9, k=2) = c * x^(-k-1) * (1-x)^(k+1)
a, b = -1e-5, 1  # domain
# note that I first change variable to the [0, 1] then do the importance sampling
d = truncated(Normal(1e-3, 1e-3), 0, 1)  # distribution(after change variable, domain=[0, 1])
nodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863]
IShalton(m, a, b, d, nodesNum)


# multivariable
g = x -> 1 / (x[1] + 1) + sqrt(x[2]) + 2 * (x[3]^2) + sqrt(2 * x[4]) + cbrt(x[5])
A, B = [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]  # domain
nodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863]
MCM(g, A, B, nodesNum, seed=1234)
quaMCM(g, A, B, nodesNum)
GLegendre(g, A, B, 30)
```


### note
- newbie of julia programming but with enthusiasm
- still in progress...
