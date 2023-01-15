# Basic Integraion
*pick the right method and don't need to worry about the domain*

## introduction
This project provide and interface for solving integration problems using guassian quadrature and monte carlo integratioon no matter the problem is univariable or multivarialbe. Furthermore, based on interest, I also implement the function for importance samplingin quasi way through halton sequences.

For Gaussian Quadrature method I implement:
1. Gauss Hermite
2. Gauss Laguerre
3. Gauss Legendre

**note: rely on this project [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) to generate the nodes and weights required by the gaussian quadrature method**

As for monte carlo method I implement:
1. naive monte carlo
2. quasi monte carlo
3. importance samplng
since the quais montecarlo method rely on the halton sequences with the prime = 2, it's better to let the number of nodes = $2^n−1$

**note: both the quasi monte carlo and importance sampling use [HaltonSequences.jl](https://github.com/tobydriscoll/HaltonSequences.jl) to implment**



## function
- `GHermite()`, `GLaguerre()`, `GLegendre()`, `GQ()`
    - `GQ()` takes advantage of gausshermite and gausslaguerre to deal with univariable problem with domain [-Inf, 0] 
    - `GQ()` will automatically apply different Gauss Quadraturer rules based on diffferent domain in multivariable problem 
        - when domain = [-Inf, Inf], use Gauss-Hermite
        - when domain = [a, Inf], use Gauss-Laguerre
        - when domain = [-Inf, b], use Gauss-legendre(it might not be the best though)
        - when the domain is proper, use Gauss-legendre
    - all the others can estimate multivariable problem with the vector input 
- `MCM()`, `quaMCM()`, `IShalton()`
    - `IShalton()` can only be used for univariable problem just now


## usage
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
