include("lzz_p.jl")

mutable struct zz_pX{T}
    _rep::vec_zz_p{T}
    function zz_pX{T}(X::Int) where {T}
        x = convert(zz_p{T}, X)
        # y= [x]
        # println("-------->", x,"----->",y)
        # println(typeof(y))
        new( [x])
    end
    function zz_pX{T}(x::zz_p{T})  where {T}
        new([x])
    end
end
