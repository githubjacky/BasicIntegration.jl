### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ cd98b09f-9bed-4255-a5f5-3b88e027a0ea
using Pkg; Pkg.activate("/Users/jacky/Documents/project/BasicIntegration.jl/")

# ╔═╡ dc9da99e-949b-11ed-24b5-c972c2336f63
begin
	using StatsPlots, HaltonSequences, Random, QuadGK, Distributions, StaticArrays
	using BasicIntegration
	g(x; c=1e-9, k=2) = c * x^(-k-1) * (1-x)^(k+1)
	
	a = 1e-5  # upper bound
	b = 1  # lower bound
	
	trueVal = quadgk(g, 1e-5, 1)[1]
	# n = 19, 20, 21, 23, 26
	NodesNum = [524_287, 1_048_575, 2_097_151, 8_388_607, 67_108_863];
end

# ╔═╡ 5095fcd3-0770-460a-88d8-29debc72d965
plot(g, xlim=[0.01, 1], label="", xlabel="x", ylabel="g(x)")

# ╔═╡ 89f9eb9d-dadc-4da0-bf91-a16cee5a3fff
GL_res = GLegendre(g, a, b, 3000)

# ╔═╡ 768c6d30-1e54-481e-99f4-cba295f0f467
MCM_deviate = abs.(MCM(g, a, b, NodesNum, seed=20220923) .- trueVal)

# ╔═╡ 24bb251f-8094-4215-841c-df038e9e4b68
quaMCM_deviate = abs.(quaMCM(g, a, b, NodesNum) .- trueVal)

# ╔═╡ 24cb3ee5-c7f0-4cf7-9e36-9ae464474018
begin
	d = truncated(Normal(1e-3, 1e-3), 0, 1)
	IS_deviate = abs.(IShalton(g, a, b, d, NodesNum) .- trueVal)
end

# ╔═╡ 91c204e1-59f8-4905-96bb-469b99cdd241
begin
	deviateMatrix = [MCM_deviate quaMCM_deviate IS_deviate]
	x = ["    0.5    ", "   1   ", "  2  ", " 8.3 ", "67.1"]
	x = repeat(x, outer=5)
	label = ["   MCM   ", "  Quasi MCM  ", "IS_halton"]
	label = repeat(label, inner=5)
	groupedbar(x, deviateMatrix, group=label, xlabel="nodes Used(million)",
	            ylabel="deviation(absolute value)", title="deviation comparison",
	            bar_width=0.5, lw=0, framestyle=:box)
end

# ╔═╡ Cell order:
# ╠═cd98b09f-9bed-4255-a5f5-3b88e027a0ea
# ╠═dc9da99e-949b-11ed-24b5-c972c2336f63
# ╠═5095fcd3-0770-460a-88d8-29debc72d965
# ╠═89f9eb9d-dadc-4da0-bf91-a16cee5a3fff
# ╠═768c6d30-1e54-481e-99f4-cba295f0f467
# ╠═24bb251f-8094-4215-841c-df038e9e4b68
# ╠═24cb3ee5-c7f0-4cf7-9e36-9ae464474018
# ╠═91c204e1-59f8-4905-96bb-469b99cdd241
