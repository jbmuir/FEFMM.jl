function mul_analytic!(τ::Array{R,1}, dx1::R, ie1::Int, xs1::Int) where {R <: Real}
    for i1 = 1:ie1
        @inbounds τ[i1] *= sqrt(sum(((i1-xs1)*dx1)^2))
    end
end

function mul_analytic!(τ::Array{R,2}, dx1::R, dx2::R, ie1::Int, ie2::Int, xs1::Int, xs2::Int) where {R <: Real}
    for i2 = 1:ie2
        for i1 = 1:ie1
            @inbounds τ[i1,i2] *= sqrt(sum(((i1-xs1)*dx1)^2+((i2-xs2)*dx2)^2))
        end
    end
end

function mul_analytic!(τ::Array{R,3}, dx1::R, dx2::R, dx3::R, ie1::Int, ie2::Int, ie3::Int, xs1::Int, xs2::Int, xs3::Int) where {R <: Real}
    for i3 = 1:ie3
        for i2 = 1:ie2
            for i1 = 1:ie1
                @inbounds τ[i1,i2,i3] *= sqrt(sum(((i1-xs1)*dx1)^2+((i2-xs2)*dx2)^2+((i3-xs3)*dx3)^2))
            end
        end
    end
end

function grad_analytic(τ0::Array{R,1}, dx1::R, ie1::Int, xs1::Int) where {R <: Real}
    ∇τ0 = [ones(size(τ0))]
    for i1 = 1:ie1
        @inbounds ∇τ0[1][i1] = dx1*(i1-xs1)/τ0[i1]
    end
    ∇τ0
end

function grad_analytic(τ0::Array{R,2}, dx1::R, dx2::R, ie1::Int, ie2::Int, xs1::Int, xs2::Int) where {R <: Real}
    ∇τ0 = [ones(size(τ0)), ones(size(τ0))]
    for i2 = 1:ie2
        for i1 = 1:ie1
            @inbounds ∇τ0[1][i1,i2] = dx1*(i1-xs1)/τ0[i1,i2]
            @inbounds ∇τ0[2][i1,i2] = dx2*(i2-xs2)/τ0[i1,i2]
        end
    end
    ∇τ0
end

function grad_analytic(τ0::Array{R,3}, dx1::R, dx2::R, dx3::R, ie1::Int, ie2::Int, ie3::Int, xs1::Int, xs2::Int, xs3::Int) where {R <: Real}
    ∇τ0 = [ones(size(τ0)), ones(size(τ0)), ones(size(τ0))]
    for i3 = 1:ie3
        for i2 = 1:ie2
            for i1 = 1:ie1
                @inbounds ∇τ0[1][i1,i2,i3] = dx1*(i1-xs1)/τ0[i1,i2,i3]
                @inbounds ∇τ0[2][i1,i2,i3] = dx2*(i2-xs2)/τ0[i1,i2,i3]
                @inbounds ∇τ0[3][i1,i2,i3] = dx3*(i3-xs3)/τ0[i1,i2,i3]
            end
        end
    end
    ∇τ0
end