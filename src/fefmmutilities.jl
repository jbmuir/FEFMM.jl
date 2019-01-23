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

function set_neighbours!(xn::Array{S,1}, x::S, tags::Array{UInt8}, cs::Array{S,1}, I1::S, Iend::S) where {S <: CartesianIndex}
    for (i,s) in enumerate(cs)
        if cilt(I1, x-s) && tags[x-s] < 0x3 
            xn[2*i-1] = x-s
        else
            xn[2*i-1] = x
        end
        if cilt(x+s, Iend) && tags[x+s] < 0x3
            xn[2*i] = x+s
        else
            xn[2*i] = x
        end
    end
end


function check1bwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, I1::S) where {R <: Real, S <: CartesianIndex}
    if cilt(I1, x-s) && tags[x-s] == 0x3
        return true
    else
        return false
    end
end

function check1fwd(τ1::Array{R}, τ0::Array{R}, tags::Array{UInt8}, x::S, s::S, Iend::S) where {R <: Real, S <: CartesianIndex}
    if cilt(x+s, Iend) && tags[x+s] == 0x3
        return true
    else
        return false
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

# function check1bwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, i::Int, Ipre::S1, Ipost::S2) where {R <: Real, N, S1 <: CartesianIndex, S2 <: CartesianIndex}
#     if 1 <= i-1 && (@inbounds tags[Ipre, i-1, Ipost] == 0x3)
#         return true
#     else
#         return false
#     end
# end

# function check2bwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, i::Int, iend::Int, Ipre::S1, Ipost::S2) where {R <: Real, N, S1 <: CartesianIndex, S2 <: CartesianIndex}
#     if 1 <= i-2 && (@inbounds tags[Ipre, i-2, Ipost] == 0x3) && (τ1[Ipre, i-1, Ipost]*τ0[Ipre, i-1, Ipost] >= τ1[Ipre, i-2, Ipost]*τ0[Ipre, i-2, Ipost])
#         return true
#     else
#         return false
#     end
# end

# function check1fwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, i::Int, iend::Int, Ipre::S1, Ipost::S2) where {R <: Real, N, S1 <: CartesianIndex, S2 <: CartesianIndex}
#     if i+1 <= iend && (@inbounds tags[Ipre, i+1, Ipost] == 0x3)
#         return true
#     else
#         return false
#     end
# end

# function check2fwd(τ1::Array{R, N}, τ0::Array{R, N}, tags::Array{UInt8, N}, i::Int, iend::Int, Ipre::S1, Ipost::S2) where {R <: Real, N, S1 <: CartesianIndex, S2 <: CartesianIndex}
#     if i+1 <= iend && (@inbounds tags[Ipre, i+2, Ipost] == 0x3) && (τ1[Ipre, i+1, Ipost]*τ0[Ipre, i+1, Ipost] > τ1[Ipre, i+2, Ipost]*τ0[Ipre, i+2, Ipost])
#         return true
#     else
#         return false
#     end
# end

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

function solve_piecewise_quadratic(α1::R, β1::R, κ2x::R) where {R <: Real}
    #solve_quadratic(α1^2, -α1*β1, β1^2-κ2x) = 
    (sqrt(κ2x)+β1)/α1
end

function solve_piecewise_quadratic(α1::R, α2::R, β1::R, β2::R, κ2x::R) where {R <: Real}
    a = α1^2+α2^2
    bo2 = -α1*β1-α2*β2
    c = β1^2+β2^2-κ2x
    τ1t = solve_quadratic(a, bo2, c)
    ∇τ1t1 = α1*τ1t-β1
    ∇τ1t2 = α2*τ1t-β2
    if (∇τ1t1 <= 0) || (∇τ1t2 <= 0)
        tr1 = abdiv(α1, β1)
        tr2 = abdiv(α2, β2)
        if tr1 <= tr2 
            τ1t = solve_piecewise_quadratic(α1, β1, κ2x)
        else
            τ1t = solve_piecewise_quadratic(α2, β2, κ2x)
        end
    end
    τ1t 
end

function solve_piecewise_quadratic(α1::R, α2::R, α3::R, β1::R, β2::R, β3::R, κ2x::R) where {R <: Real}
    a = α1^2+α2^2+α3^2
    bo2 = -α1*β1-α2*β2-α3*β3
    c = β1^2+β2^2+β3^2-κ2x
    τ1t = solve_quadratic(a, bo2, c)
    ∇τ1t1 = α1*τ1t-β1
    ∇τ1t2 = α2*τ1t-β2
    ∇τ1t3 = α3*τ1t-β3
    if (∇τ1t1 <= 0) || (∇τ1t2 <= 0) || (∇τ1t3 <= 0)
        tr1 = abdiv(α1, β1)
        tr2 = abdiv(α2, β2)
        tr3 = abdiv(α3, β3)
        if tr1 <= tr2 
            if tr2 <= tr3
                τ1t = solve_piecewise_quadratic(α1, α2, β1, β2, κ2x)
            else
                τ1t = solve_piecewise_quadratic(α1, α3, β1, β3, κ2x)
            end
        else
            if tr1 <= tr3
                τ1t = solve_piecewise_quadratic(α2, α1, β2, β2, κ2x)
            else
                τ1t = solve_piecewise_quadratic(α2, α3, β2, β3, κ2x)
            end
        end
    end
    τ1t 
end

# function solve_piecewise_quadratic(α::Array{R}, β::Array{R}, κ2x::R) where {R <: Real}
#     a = sum(α.^2)
#     bo2 = -sum((α.^2).*β)
#     c = sum((α.^2).*(β.^2))-κ2x
#     τ1t = solve_quadratic(a, bo2, c)
#     ∇τ1t = α.*(τ1t.-β) 
#     #this loop gets rid of the value with the highest beta / alpha - ie that with the arithmetically smallest estimated derivative
#     if any(∇τ1t .<= zero(R))
#         tr = abdiv.(α,β)
#         p = sortperm(tr)[1:(end-1)]
#         τ1t = solve_piecewise_quadratic(α[p], β[p], κ2x)
#     end
#     τ1t  
# end
