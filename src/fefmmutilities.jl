struct Node{R <: Real}
    ind::CartesianIndex
    val::R
end

<(a::Node, b::Node) = a.val < b.val

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

function check1bwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, I1::S) where {R <: Real, S <: CartesianIndex}
    if cilt(I1, x-s) && tags[x-s] == 0x3
        return (true, τ1[x-s]*τ0[x-s])
    else
        return (false, -Inf)
    end
end

function check1fwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, Iend::S) where {R <: Real, S <: CartesianIndex}
    if cilt(x+s, Iend) && tags[x+s] == 0x3
        return (true, τ1[x+s]*τ0[x+s])
    else
        return (false, -Inf)
    end
end

function check2bwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, I1::S) where {R <: Real, S <:CartesianIndex}
    if cilt(I1, x-2*s) && tags[x-2*s] == 0x3 && τ1[x-s]*τ0[x-s] >= τ1[x-2*s]*τ0[x-2*s]
        return true
    else
        return false
    end 
end

function check2fwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, Iend::S) where {R <: Real, S <:CartesianIndex}
    if cilt(x+2*s, Iend) && tags[x+2s] == 0x3 && τ1[x+s]*τ0[x+s] > τ1[x+2*s]*τ0[x+2*s]
        return true
    else
        return false
    end 
end

function neighbours(x::S, tags::Array{UInt8}, cs::Array{S}, I1::S, Iend::S) where {S <: CartesianIndex}
    xn = Array{S}(undef, 0)
    for s in cs
        if cilt(I1, x-s) && tags[x-s] < 0x3 
            push!(xn, x-s)
        end
        if cilt(x+s, Iend) && tags[x+s] < 0x3
            push!(xn, x+s)
        end
    end
    xn
end

function mul_analytic!(τ0::Array{R}, dx::Array{R}, x0::CartesianIndex)  where {R <: Real}
    inds = CartesianIndices(τ0)
    for i in inds
        τ0[i] *= sqrt(sum(((Tuple(i).-Tuple(x0)).*dx).^2))
    end
end

function mul_grad_analytic!(τ0::Array{R}, dx::Array{R}, x0::CartesianIndex, dim::Int)  where {R <: Real}
    inds = CartesianIndices(τ0)
    for i in inds
        τ0[i] *= dx[dim].*(Tuple(i)[dim]-Tuple(x0)[dim])./sqrt(sum(((Tuple(i).-Tuple(x0)).*dx).^2))
    end
end

function solve_quadratic(a::R, bo2::R, c::R) where {R <: Real}
    #bo2 = b over 2
    in_sqrt = bo2^2 - a*c
    if in_sqrt > zero(R)
        return (-bo2 + sqrt(in_sqrt))/a
    else
        return zero(R)
    end
end

function abdiv(a::R, b::R) where {R <: Real}
    if a == zero(R)
        return convert(R, Inf)
    else
        return b/a
    end
end

function solve_piecewise_quadratic(α::Array{R}, β::Array{R}, κ2x::R) where {R <: Real}
    a = sum(α.^2)
    bo2 = -sum((α.^2).*β)
    c = sum((α.^2).*(β.^2))-κ2x
    τ1t = solve_quadratic(a, bo2, c)
    ∇τ1t = α.*(τ1t.-β) 
    #this loop gets rid of the value with the highest beta / alpha - ie that with the arithmetically smallest estimated derivative
    if any(∇τ1t .<= zero(R))
        tr = abdiv.(α,β)
        p = sortperm(tr)[1:(end-1)]
        τ1t = solve_piecewise_quadratic(α[p], β[p], κ2x)
    end
    τ1t  
end