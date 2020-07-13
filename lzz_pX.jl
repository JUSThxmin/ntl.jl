include("lzz_p.jl")
import Base: getindex, zero, one


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
    function zz_pX{T}() where {T}
        return zz_pX{T}(0)
    end
    function zz_pX{T}(x::zz_p{T})  where {T}
        new([x])
    end
    function zz_pX{T}(v::vec_zz_p{T}) where {T}
        p=new(v)
        normalize!(p)
    end
    function zz_pX{T}(v::Array{Int64,1}) where {T}
        p=new(v)
        normalize!(p)
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
function normalize!(p::zz_pX{T})::zz_pX{T} where {T}
    r = p._rep
    m= length(r)
    while (m > 1) & iszero(r[m])
        pop!(r)
        m -= 1
    end
    return p
end

function deg(p::zz_pX{T}) where {T}
    return length(p._rep) - 1
end
function coeff(p::zz_pX{T}, i::Int) where {T}
    if i < 0
        error("coeff: negative index")
    end
    if i in 0:deg(p)
        return p._rep[i+1]
    end
    return zz_p{T}(0)
end

function zero(::Type{zz_pX{T}}) where {T}
    return zz_p{T}(0)
end
function one(::Type{zz_pX{T}}) where {T}
    return zz_p{T}(1)
end

function GetCoeff!(x::zz_p{T}, p::zz_pX{T}, i::Int) where {T}
    x = coeff(p,i)
end
function LeadCoeff(p::zz_pX{T}) where {T} 
    p._rep[end] 
end
function ConstTerm(p::zz_pX{T}) where {T}
    p._rep[1]
end
function getindex(p::zz_pX{T}, i::Int) where {T}
    if i < 0
        error("getindex: negative index")
    end
    if i in 0:deg(p)
        return getindex(p._rep, i+1)
    end
    return zero(zz_pX{T})
end

function setindex!(p::zz_pX{T}, i::Int, x::zz_p{T}) where {T}
    if i < 0
        error("setindex!: nagative index")
    end
    m = deg(p)
    if (i > m) & !(iszero(x)) 
        return x
    end
    if i >m
        z=zeros(zz_p{T},i-m)
        append!(p._rep, z)
    end
    setindex!(p._rep, i, x)
    normalize!(p)
end


