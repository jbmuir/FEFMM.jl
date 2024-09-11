# function schecklt(a::Int, b1::Int, b2::Int, s::Int)
#     a < b && cld(b1, s) == cld(b2, s)
# end

cilt(a::S, b::S) where {S <: CartesianIndex} = all(a.I .<= b.I)

function cartstrides(A::AbstractArray)
    s = size(A)
    dims = length(s)
    inds = []
    for i = 1:dims
        z = zeros(Int, dims)
        z[i] = 1
        push!(inds, CartesianIndex(z...))
    end
    [inds...]
end

function set_neighbours!(xn::Vector{S}, x::S, tags::Array{UInt8, N}, cs::Vector{S}, I1::S, Iend::S) where {N, S <: CartesianIndex{N}}
    for (i,s) in enumerate(cs)
        if cilt(I1, x-s) && @inbounds tags[x-s] < 0x3 
            @inbounds xn[2*i-1] = x-s
        else
            @inbounds xn[2*i-1] = x
        end
        if cilt(x+s, Iend) && @inbounds tags[x+s] < 0x3
            @inbounds xn[2*i] = x+s
        else
            @inbounds xn[2*i] = x
        end
    end
end


@noinline function check1bwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, x::S, s::S, I1::S) where {N, R <: AbstractFloat, S <: CartesianIndex{N}}
    xms = x - s
    return cilt(I1, xms) && @inbounds tags[xms] == 0x3
end

@noinline function check1fwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, x::S, s::S, Iend::S) where {N, R <: AbstractFloat, S <: CartesianIndex{N}}
    xps = x + s
    return cilt(xps, Iend) && @inbounds tags[xps] == 0x3
end

@noinline function check2bwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, x::S, s::S, I1::S) where {N, R <: AbstractFloat, S <:CartesianIndex{N}}
    xms = x - s
    xm2s = x - 2*s
    return cilt(I1, xm2s) && @inbounds tags[xm2s] == 0x3 && @inbounds τ1[xms]*τ0[xms] >= τ1[xm2s]*τ0[xm2s]
end

@noinline function check2fwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, x::S, s::S, Iend::S) where {N, R <: AbstractFloat, S <:CartesianIndex{N}}
    xps = x + s
    xp2s = x + 2*s
    return cilt(xp2s, Iend) && @inbounds tags[xp2s] == 0x3 && @inbounds τ1[xps]*τ0[xps] > τ1[xp2s]*τ0[xp2s]
end

