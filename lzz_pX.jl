include("lzz_p.jl")
import Base: setindex!, getindex, zero, one, copy,>>, <<


## vector of zz_p{T}
const vec_zz_p{T} = Vector{zz_p{T}} where {T} 

function convert(::Type{vec_zz_p{T}}, v::Vector{Int}) where {T}
    return map(zz_p{T}, v)
end

## polynomial with coefficients in zz_p{T}
struct NTL_INIT_zz_pX end
ntl_init_zz_pX = NTL_INIT_zz_pX()

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
    function zz_pX{T}(v::Array{Int64,1},::NTL_INIT_zz_pX) where {T}
        new(v)
    end
    function zz_pX{T}(v::vec_zz_p{T},::NTL_INIT_zz_pX) where {T}
        new(v)
    end
end


convert(::Type{zz_pX{T}}, X::Int) where {T} = zz_pX{T}(X)
convert(::Type{zz_pX{T}}, X::zz_p{T}) where {T} = zz_pX{T}(X)

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
    while (m > 1) && iszero(r[m])
        pop!(r)
        m -= 1
    end
    return p
end

@inline function deg(p::zz_pX{T}) where {T}
    return length(p._rep)-1
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
    return zz_pX{T}(0)
end
function one(::Type{zz_pX{T}}) where {T}
    return zz_pX{T}(1)
end
function setX(::Type{zz_pX{T}}) where {T}
    zz_pX{T}([0,1],ntl_init_zz_pX)
end

@inline function isX(x::zz_pX{T}) where {T}
   return (deg(x) == 1) && isone(x._rep[2]) && iszero(x._rep[1]);
end
@inline function iszero(x::zz_pX{T}) where {T}
    return (deg(x)== 0) && iszero(x._rep[1]) 
end
@inline function isone(x::zz_pX{T}) where {T}
    return (deg(x)== 0) && isone(x._rep[1]) 
end

function GetCoeff!(x::zz_p{T}, p::zz_pX{T}, i::Int) where {T}
    x._rep = coeff(p,i)._rep
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

function setindex!(p::zz_pX{T},  x::zz_p{T}, i::Int) where {T}
    if i < 0
        error("setindex!: nagative index")
    end
    m = deg(p)
    if (i > m) && (iszero(x)) 
        return x
    end
    if i>m
        z=zeros(zz_p{T},i-m)
        append!(p._rep, z)
    end 
    setindex!(p._rep, x, i+1)
    normalize!(p)
end

setindex!(p::zz_pX{T},  x::Int, i::Int) where {T} = setindex!(p,zz_p{T}(x),i)

@inline function swap!(x::zz_pX{T}, y::zz_pX{T}) where {T}
    z = x._rep; x._rep = y._rep; y._rep = z
    nothing
end

"""
rand(::Type, d::Int)

generate a random polynomial of degree <= d 

## Example
julia> P7=zz_pX{7}
julia> rand(P7, 3)

"""
function rand(::Type{zz_pX{T}}, d::Int) where {T}
    if d < 0
        error("rand: nagative degree")
    end
    zz_pX{T}(rand(Int, d+1))
end

function trunc(p::zz_pX{T}, d::Int) where {T}
    if d > deg(p)
        return p
    else
        zz_pX{T}(p._rep[1:d])
    end
end
function trunc!(p::zz_pX{T}, d::Int) where {T}
    if d > deg(p)
        return p
    else
        resize!(p._rep,d)
    end
    normalize!(p)
end
"""
RightShift!(p::zz_pX{T}, d::Int)
p /= x^d

"""
@inline function RightShift!(p::zz_pX{T}, d::Int) where {T}
    if d>deg(p)
        p._rep = [ zero(zz_p{T}) ]
    else
        deleteat!(p._rep, 1:d)
    end
    return p
end
@inline function RightShift(p::zz_pX{T}, d::Int) where {T}
    if d>deg(p)
        return zero(zz_pX{T})
    else
        return zz_pX{T}(p._rep[1+d:end],ntl_init_zz_pX)
    end
 end
 >>(p::zz_pX{T}, d::Int) where {T} = RightShift(p,d)

"""
LeftShift!(p::zz_pX{T}, d::Int)
p *= x^d

"""
@inline function LeftShift!(p::zz_pX{T}, d::Int) where {T}
    if (d>0) && !(iszero(p))
        prepend!(p._rep,zeros(zz_p{T}, d))
    end
    return p
end
@inline function LeftShift(p::zz_pX{T}, d::Int) where {T}
    LeftShift!(copy(p),d)
end
<<(p::zz_pX{T}, d::Int) where {T} = LeftShift(p,d)

function copy(p::zz_pX{T}) where {T}
    zz_pX{T}(copy(p._rep), ntl_init_zz_pX)
end

==(a::zz_pX{T}, b::zz_pX{T}) where {T} = a._rep == b._rep

@inline function diff(a::zz_pX{T}) where {T} 
    b=a>>1
    b._rep.*=collect(1:length(b._rep))
    normalize!(b)
end


# ***************************************************************
#                          Addition
# ***************************************************************
@inline function add(x::zz_pX{T}, y::zz_pX{T}) where {T}
    dx =length(x._rep); dy = length(y._rep)
    if dx>dy
        _add_large(x,dx,y, dy)
    elseif dx<dy
        _add_large(y,dy,x, dx)
    else
        _add_equal(x,y)
    end
end
@inline function _add_equal(x::zz_pX{T}, y::zz_pX{T}) where{T}
    zz_pX{T}(x._rep + y._rep)
end
@inline function _add_large(x::zz_pX{T} , dx::Int, y::zz_pX{T}, dy::Int ) where {T}
    z=append!(x._rep[1:dy]+y._rep, x._rep[dy+1:dx])
    zz_pX{T}(z, ntl_init_zz_pX)
end
+(x::zz_pX{T},y::zz_pX{T}) where {T} = add(x, y)
add(X::zz_p{T},y::zz_pX{T}) where {T} = add(convert(zz_pX{T},X), y)
add(X::Int,y::zz_pX{T}) where {T} = add(convert(zz_pX{T},X), y)
add(x::zz_pX{T},Y::zz_p{T}) where {T} = add(x,convert(zz_pX{T},Y))
add(x::zz_pX{T},Y::Int) where {T} = add(x,convert(zz_pX{T},Y))

