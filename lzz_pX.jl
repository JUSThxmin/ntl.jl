include("lzz_p.jl")

## vector of zz_p{T}
vec_zz_p{T} = Vector{zz_p{T}} where {T} 

function convert(::Type{vec_zz_p{T}}, v::Vector{Int}) where {T}
    return map(zz_p{T}, v)
end

## polynomial with coefficients in zz_p{T}

mutable struct zz_pX{T}
    _rep::vec_zz_p{T}
    function zz_pX{T}(X::Int) where {T}
        x = zz_p{T}(X)
        new([x])
    end
    function zz_pX{T}(x::zz_p{T})  where {T}
        new([x])
    end
    function zz_pX{T}(v::vec_zz_p{T}) where {T}
        new(v)
    end
    function zz_pX{T}(v::Array{Int64,1}) where {T}
        new(v)
    end
end

"""
show(io::IO, p::zz_pX{T}, x::Symbol) 

show the polynomial p in variable x
julia> p = zz_pX{7}([1,2,3])
julia> show(stdout, p, :x)
"""
function show(io::IO, p::zz_pX{T}, sym::Symbol) where {T}
    coef = p._rep
    print(io, coef[1])
    x = string(sym)
    len = length(p._rep)
    for i = 2:len
        if !iszero(coef[i])
        print(io, "+", coef[i], "*", x, "^", i-1)
        end
    end
    return nothing
end
