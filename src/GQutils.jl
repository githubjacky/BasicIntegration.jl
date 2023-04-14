function GHdt(a, b)
    if b == Inf
        x1, jcb1 = GQrule2(a)
    elseif a == -Inf
        x1, jcb1 = GQrule3(b)
    elseif a != -1 || b != 1
        x1, jcb1 = GQrule4(a, b)
    else
        x1, jcb1 = t -> t, t -> 1
    end
    x2, jcb2 = GQrule5()
    return x1, jcb1, x2, jcb2
end


function MultiDim_GHdt(g::Function, A::Vector{T}, B::Vector{P}, invp::Function, dim) where {T<:Real, P<:Real}
    function trans(t) 
        Jcb = 1
        for i = 1:dim
            a, b = A[i], B[i]
            if a == -Inf && b == Inf
                x1, jcb1, x2, jcb2 = t->t, t->1, t->t, t->1
            else
                x1, jcb1, x2, jcb2 = GHdt(a, b)
            end
            Jcb *= (invp(t[i]) * jcb1(x2(t[i])) * jcb2(t[i]))
            t[i] = x1(x2(t[i]))
        end
        return g(t) * Jcb
    end
    return trans
end


function GLadt(a, b)
    if b == Inf
        x1, jcb1 = t -> t + a, t -> 1
        x2, jcb2 = t -> t, t -> 1
        return x1, jcb1, x2, jcb2
    elseif a == -Inf && b == Inf
        x1, jcb1 = GQrule1()
    elseif a == -Inf
        x1, jcb1 = GQrule3(b)
    elseif a != -1 || b != 1
        x1, jcb1 = GQrule4(a, b)
    else
        x1, jcb1 = t -> t, t -> 1
    end
    x2, jcb2 = GQrule6()
    return x1, jcb1, x2, jcb2
end

function MultiDim_GLadt(g::Function, A::Vector{T}, B::Vector{P}, invp::Function, dim) where {T<:Real, P<:Real}
    function trans(t) 
        Jcb = 1
        for i = 1:dim
            a, b = A[i], B[i]
            if a == 0 && b == Inf
                x1, jcb1, x2, jcb2 = t->t, t->1, t->t, t->1
            else
                x1, jcb1, x2, jcb2 = GLadt(a, b)
            end
            Jcb *= (invp(t[i]) * jcb1(x2(t[i])) * jcb2(t[i]))
            t[i] = x1(x2(t[i]))
        end
        return g(t) * Jcb
    end
    return trans
end


function GLdt(a, b)
    if a == -Inf && b == Inf
        x1, jcb1 = GQrule1()
    elseif b == Inf
        x1, jcb1 = GQrule2(a)
    elseif a == -Inf
        x1, jcb1 = GQrule3(b)
    else
        x1, jcb1 = GQrule4(a, b)
    end
    return x1, jcb1
end

function MultiDim_GLdt(g::Function, A::Vector{T}, B::Vector{P}, dim) where {T<:Real, P<:Real}
    function trans(t) 
        Jcb = 1
        for i = 1:dim
            a, b = A[i], B[i]
            if a == -1 && b == 1
                x1, jcb1 = t->t, t->1
            else
                x1, jcb1 = GLdt(a, b)
            end
            Jcb *= jcb1(t[i])
            t[i] = x1(t[i])
        end
        return g(t) * Jcb
    end
    return trans
end


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

function addInvp(g::Function, A::Vector{T}, B::Vector{P}, dim, groupIndex::Vector{Int64}) where {T<:Real, P<:Real}
    function trans(t)
        Invp = 1
        for i = 1:dim
            a, b = A[i], B[i]
            if a == -Inf && b == Inf
                Invp *= exp(t^2)
                groupIndex[i] = 1
            elseif b == Inf
                Invp *= exp(t)
                a != 0 && (t[i] = t[i] + a)
                groupIndex[ij = 2]
            else
                groupIndex[i] = 3
            end
        end
        return g(t) * Invp
    end
    return trans
end


function nodesGenGQ(n, dim, nodesArr::Vector{Vector{Float64}}, wArr::Vector{Vector{Float64}}, groupIndex::Vector{Int64})
    N = n^dim
    Nodes, W = Vector{Vector{Float64}}(undef, N), Vector{Float64}(undef, N)
    Nodesᵢ, Wᵢ = Vector{Float64}(undef, dim), Vector{Float64}(undef, dim)
    for i = 1:dim
        Nodesᵢ[i] = nodesArr[groupIndex[i]][1]
        Wᵢ[i] = wArr[groupIndex[i]][1]
    end
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
            Nodesᵢ[1], Wᵢ[1] = nodesArr[groupIndex[1]][index[1]], wArr[groupIndex[1]][index[1]]
        else
            for j = 2:dim
                if index[j] == n
                    index[j] = 1
                    Nodesᵢ[j], Wᵢ[j] = nodesArr[groupIndex[j]][1], wArr[groupIndex[j]][1]
                else
                    index[j] += 1
                    Nodesᵢ[j], Wᵢ[j] = nodesArr[groupIndex[j]][index[j]], wArr[groupIndex[j]][index[j]]
                    break
                end
            end
            index[1] = 1
            Nodesᵢ[1], Wᵢ[1] = nodesArr[groupIndex[1]][1], wArr[groupIndex[1]][1]
        end
    end
    return Nodes, W
end


function GQrule1()  # [-1, -1] to [-Inf, Inf]
    x(t) = t / (1 - t^2)
    jcb(t) = (t^2 + 1) / (t^2 - 1)^2
    return x, jcb
end

function GQrule2(a)  # [-1, 1] to [a, Inf]
    x(t) = a + (1 + t) / (1 - t)
    jcb(t) = 2 / (t - 1)^2
    return x, jcb
end

function GQrule3(b) # [-1, 1] to [-Inf, b]
    x(t) = b - (1 - t) / (1 + t)
    jcb(t) = 2 / (t + 1)^2
    return x, jcb
end

function GQrule4(a, b) # [-1, 1] to proper
    x(t) = (a + b) / 2 + (b - a) * t / 2
    jcb(t) = (b - a) / 2
    return x, jcb
end

function GQrule5()  # [-Inf, Inf] to [-1, 1]
    x(t) = (-1 + sqrt(1 + 4 * t^2)) / (2 * t)
    jcb(t) = ((8 * t^2 * (sqrt(1 + 4 * t^2))^(-1)) -
              (2 * (-1 + sqrt(1 + 4 * t^2)))) / (4 * t^2)
    return x, jcb
end

function GQrule6()  # [0, Inf] to [-1, 1]
    x(t) = (t - 1) / (t + 1)
    jcb(t) = (2 / (t + 1)^2)
    return x, jcb
end
