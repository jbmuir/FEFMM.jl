function mul_analytic!(τ::Vector{<:AbstractFloat},
                       dx1::AbstractFloat,
                       ie1::Integer,
                       xs1::Integer)
    for i1 = 1:ie1
        @inbounds τ[i1] *= sqrt(sum(((i1-xs1)*dx1)^2))
    end
end

function mul_analytic!(τ::Matrix{<:AbstractFloat},
                       dx1::AbstractFloat, dx2::AbstractFloat,
                       ie1::Integer, ie2::Integer,
                       xs1::Integer, xs2::Integer)
    for i2 = 1:ie2
        for i1 = 1:ie1
            @inbounds τ[i1,i2] *= sqrt(sum(((i1-xs1)*dx1)^2+((i2-xs2)*dx2)^2))
        end
    end
end

function mul_analytic!(τ::Array{<:AbstractFloat,3},
                       dx1::AbstractFloat, dx2::AbstractFloat, dx3::AbstractFloat,
                       ie1::Integer, ie2::Integer, ie3::Integer,
                       xs1::Integer, xs2::Integer, xs3::Integer)
    for i3 = 1:ie3
        for i2 = 1:ie2
            for i1 = 1:ie1
                @inbounds τ[i1,i2,i3] *= sqrt(sum(((i1-xs1)*dx1)^2+((i2-xs2)*dx2)^2+((i3-xs3)*dx3)^2))
            end
        end
    end
end

function grad_analytic(τ0::Vector{R},
                       dx1::AbstractFloat,
                       ie1::Integer,
                       xs1::Integer) where {R <: AbstractFloat}
    ∇τ0 = [ones(R, size(τ0))]
    for i1 = 1:ie1
        @inbounds ∇τ0[1][i1] = dx1*(i1-xs1)/τ0[i1]
    end
    #set terms at source
    ∇τ0[1][xs1] = 1
    ∇τ0
end

function grad_analytic(τ0::Matrix{R},
                       dx1::AbstractFloat, dx2::AbstractFloat,
                       ie1::Integer, ie2::Integer,
                       xs1::Integer, xs2::Integer) where {R <: AbstractFloat}
    ∇τ0 = [ones(R, size(τ0)), ones(R, size(τ0))]
    for i2 = 1:ie2
        for i1 = 1:ie1
            @inbounds ∇τ0[1][i1,i2] = dx1*(i1-xs1)/τ0[i1,i2]
            @inbounds ∇τ0[2][i1,i2] = dx2*(i2-xs2)/τ0[i1,i2]
        end
    end
    ∇τ0[1][xs1, xs2] = 1/sqrt(2)
    ∇τ0[2][xs1, xs2] = 1/sqrt(2)
    ∇τ0
end

function grad_analytic(τ0::Array{R,3},
                       dx1::AbstractFloat, dx2::AbstractFloat, dx3::AbstractFloat,
                       ie1::Integer, ie2::Integer, ie3::Integer,
                       xs1::Integer, xs2::Integer, xs3::Integer) where {R <: AbstractFloat}
    ∇τ0 = [ones(R, size(τ0)), ones(R, size(τ0)), ones(R, size(τ0))]
    for i3 = 1:ie3
        for i2 = 1:ie2
            for i1 = 1:ie1
                @inbounds ∇τ0[1][i1,i2,i3] = dx1*(i1-xs1)/τ0[i1,i2,i3]
                @inbounds ∇τ0[2][i1,i2,i3] = dx2*(i2-xs2)/τ0[i1,i2,i3]
                @inbounds ∇τ0[3][i1,i2,i3] = dx3*(i3-xs3)/τ0[i1,i2,i3]
            end
        end
    end
    ∇τ0[1][xs1, xs2, xs3] = 1/sqrt(3)
    ∇τ0[2][xs1, xs2, xs3] = 1/sqrt(3)
    ∇τ0[3][xs1, xs2, xs3] = 1/sqrt(3)
    ∇τ0
end

