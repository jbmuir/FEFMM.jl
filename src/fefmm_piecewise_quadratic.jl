function abdiv(a::R, b::R) where {R <: AbstractFloat}
    if a == zero(R)
        return convert(R, Inf)
    else
        return b/a
    end
end

function solve_quadratic(a::R, bo2::R, c::R) where {R <: AbstractFloat}
    #bo2 = b over 2
    in_sqrt = bo2^2 - a*c
    if in_sqrt > zero(R)
        return (-bo2 + sqrt(in_sqrt))/a
    else
        return zero(R)
    end
end

function solve_piecewise_quadratic(α1::R, β1::R, κ2x::R) where {R <: AbstractFloat}
    #solve_quadratic(α1^2, -α1*β1, β1^2-κ2x) = 
    (sqrt(κ2x)+β1)/α1
end

function solve_piecewise_quadratic(α1::R, α2::R, β1::R, β2::R, κ2x::R) where {R <: AbstractFloat}
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

function solve_piecewise_quadratic(α1::R, α2::R, α3::R, β1::R, β2::R, β3::R, κ2x::R) where {R <: AbstractFloat}
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
                τ1t = solve_piecewise_quadratic(α2, α1, β2, β1, κ2x)
            else
                τ1t = solve_piecewise_quadratic(α2, α3, β2, β3, κ2x)
            end
        end
    end
    τ1t 
end

