using FastGaussQuadrature

n = 30
xi1, wi1 = gausshermite(n)


function nodesGen(n, dim, nodes::Vector{Float64}, w::Vector{Float64})
    N = n^dim
    Nodes, W = Vector{Vector{Float64}}(undef, N), Vector{Float64}(undef, N)
    Nodesᵢ, Wᵢ = repeat([nodes[1]], dim), repeat([w[1]], dim)
    index, i = ones(Int64, dim), 1
    while true
        Nodes[i] = copy(Nodesᵢ)
        u = 1
        for k in Wᵢ
            u *= k
        end
        W[i] = u
        i != N ? (i += 1) : break 
        if index[1] != n
            index[1] += 1
            Nodesᵢ[1], Wᵢ[1] = nodes[index[1]], w[index[1]]
        else
            for j = 2:dim
                if index[j] == n
                    index[j] = 1
                    Nodesᵢ[j], Wᵢ[j] = nodes[1], w[1]
                else
                    index[j] += 1
                    Nodesᵢ[j], Wᵢ[j] = nodes[index[j]], w[index[j]]
                    break
                end
            end
            index[1] = 1
            Nodesᵢ[1], Wᵢ[1] = nodes[1], w[1]
        end
    end
    return Nodes, W
end

xi2, wi2 = nodesGen(n, 1, xi1, wi1)

for i = 1:n
    xi1[i] != xi2[i][1] && (println("$i is difference"); break)
    i == n && println("all is the same")
end
