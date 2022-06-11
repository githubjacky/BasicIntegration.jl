function MCMdt(a, b)
    if a == -Inf && b == Inf
        x1, jcb1 = MCMrule1()
    elseif b == Inf
        x1, jcb1 = MCMrule2(a)
    elseif a == -Inf
        x1, jcb1 = MCMrule3(b)
    elseif a != 0 || b != 1
        x1, jcb1 = MCMrule4(a, b)
    else
        x1, jcb1 = t -> t, t -> 1
    end
    return x1, jcb1
end


function multiDim_MCMdt(g::Function, A::Vector{T}, B::Vector{P}, dim) where {T<:Real,P<:Real}
    function trans(t)
        Jcb = 1
        for i = 1:dim
            x1, jcb1 = MCMdt(A[i], B[i])
            Jcb *= jcb1(t[i])
            t[i] = x1(t[i])
        end
        return g(t) * Jcb
    end
    return trans
end


function MCMrule1()  # [-Inf, Inf] to [0, 1]
    x(t) = (2t - 1) / (t - t^2)
    jcb(t) = (2t^2 - 2t + 1) / (t^2 - t)^2
    return x, jcb
end

function MCMrule2(a)  # [a, Inf] to [0, 1]
    x(t) = a + t / (1 - t)
    jcb(t) = 1 / (t - 1)^2
    return x, jcb
end

function MCMrule3(b)  # [-Inf, b] to [0, 1]
    x(t) = b + (t - 1) / t
    jcb(t) = 1 / t^2
    return x, jcb
end

function MCMrule4(a, b)  # proper to [0, 1]
    x(t) = a + (b - a)t
    jcb(t) = b - a
    return x, jcb
end
