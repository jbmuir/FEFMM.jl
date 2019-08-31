function solve_node(τ1::Array{R, N}, 
                    α::Vector{R},
                    β::Vector{R},
                    τ0::Array{R, N},
                    ∇τ0::Array{T},
                    tags::Array{UInt8, N},
                    κ2::Array{R, N},
                    x::S,
                    dx::Vector{<:AbstractFloat},
                    cs::Vector{S},
                    I1::S, 
                    Iend::S) where {N, R <: AbstractFloat, T <: Array{R, N}, S <: CartesianIndex{N}}
    #loop over dimensions to set \alpha and \beta ...
    for (i, s) in enumerate(cs)
        fwdok = check1fwd(τ1, τ0, tags, x, s, Iend)
        bwdok = check1bwd(τ1, τ0, tags, x, s, I1)
        if fwdok && bwdok 
            if @inbounds τ1[x+s]*τ0[x+s] < τ1[x-s]*τ0[x-s]
                bwkok = false
            else 
                fwdok = false
            end
        end
        # Note Here: working out the maths gives an additional 1/\alpha for the \beta terms
        # to reduce total operations we redefine \beta = \beta' / \alpha where \beta' is the \beta in the paper
        if fwdok 
            if check2fwd(τ1, τ0, tags, x, s, Iend)
                @inbounds α[i] = 1.5*τ0[x]/dx[i]-∇τ0[i][x]
                @inbounds β[i] = (2*τ1[x+s]-0.5*τ1[x+2*s])*τ0[x]/dx[i]
            else
                @inbounds α[i] = τ0[x]/dx[i]-∇τ0[i][x]
                @inbounds β[i] = τ0[x]*τ1[x+s]/dx[i]
            end
        elseif bwdok
            if check2bwd(τ1, τ0, tags, x, s, I1)
                @inbounds α[i] = 1.5*τ0[x]/dx[i]+∇τ0[i][x]
                @inbounds β[i] = (2*τ1[x-s]-0.5*τ1[x-2*s])*τ0[x]/dx[i]
            else
                @inbounds α[i] = τ0[x]/dx[i]+∇τ0[i][x]
                @inbounds β[i] = τ0[x]*τ1[x-s]/dx[i]
            end
        else
            @inbounds α[i] = zero(R)
            @inbounds β[i] = zero(R)
        end
    end
    #All of the \alpha and \beta should be set, so now we will proceed with solving the solve_quadratic
    #@inbounds solve_piecewise_quadratic(α..., β..., κ2[x])
    if N == 1
        @inbounds τ1y = solve_piecewise_quadratic(α[1], β[1], κ2[x])
    elseif N == 2
        @inbounds τ1y = solve_piecewise_quadratic(α[1], α[2], β[1], β[2], κ2[x])
    elseif N == 3
        @inbounds τ1y = solve_piecewise_quadratic(α[1], α[2], α[3], β[1], β[2], β[3], κ2[x])        
    end
    τ1y
end


function fefmm_loop!(τ1::Array{R, N}, 
                     ocount::Integer,
                     ordering::Vector{<:Integer},
                     α::Vector{R},
                     β::Vector{R},
                     τ0::Array{R, N}, 
                     ∇τ0::Array{T},
                     tags::Array{UInt8, N}, 
                     front::BinaryMinHeap{Node{R, N}},
                     κ2::Array{R, N},
                     dx::Vector{<:AbstractFloat},
                     cs::Vector{S},
                     xn::Vector{S},
                     I1::S, 
                     Iend::S, 
                     LI::LinearIndices{N}) where {N, R <: AbstractFloat, T <: Array{R, N}, S <: CartesianIndex{N}}
    while isempty(front) == false
        x = pop!(front).ind
        @inbounds if tags[x] < 0x3# as we are using a non-mutable minheap, we might have already set this node to known on another pass (i.e. this could be an old worthless version of the node). This check makes sure we don't accidentally override something we already have computed
            @inbounds ordering[ocount] = LI[x]
            ocount += 1
            @inbounds tags[x] = 0x3
            set_neighbours!(xn, x, tags, cs, I1, Iend)
            for y in xn 
                if y != x # to avoid allocating neighbours repeatedly we reuse xn but have y = x as the "neighbour not acceptable" flag value
                    @inbounds tags[y] = 0x2 #move neighbours to front
                    τ1y = solve_node(τ1, α, β, τ0, ∇τ0, tags, κ2, y, dx, cs, I1, Iend)
                    @inbounds if τ1y < τ1[y]
                            #only update and add new node to heap if the new solution is smaller 
                        @inbounds τ1[y] = τ1y
                        @inbounds push!(front, Node(y, τ1[y]*τ0[y]))
                    end
                end
            end
        end
    end
end


"""
Inputs to fefmm: 
κ - array of slowness
dx - array of grid spacings
x0 - linear index of source location

Outputs of fefmm: 
τ1τ0 = travel time array from source

Tags: 
1 = unknown
2 = front
3 = known
"""

function fefmm(κ2::Array{R, N},
               dx::Vector{<:AbstractFloat},
               xs::CartesianIndex{N}) where {R <: AbstractFloat, N}
    #initialization
    ordering = Array{Int}(undef, length(κ2))
    cs = cartstrides(κ2)
    inds = CartesianIndices(κ2)
    I1 = first(inds)
    Iend = last(inds)
    τ0 = ones(R, size(κ2))
    mul_analytic!(τ0, dx..., size(κ2)..., Tuple(xs)...)
    ∇τ0 = grad_analytic(τ0, dx..., size(κ2)..., Tuple(xs)...)
    τ1 = R(Inf).*ones(R, size(κ2))
    τ1[xs] = sqrt(κ2[xs])
    tags = ones(UInt8,size(κ2))
    tags[xs] = 0x2
    front = BinaryMinHeap{Node{R,N}}()
    push!(front, Node(inds[xs], τ1[xs]*τ0[xs])) 
    α = Array{R}(undef,length(size(τ1)))
    β = similar(α)
    LI = LinearIndices(τ1)
    xn = Array{typeof(I1)}(undef, N*2)
    #main loop
    fefmm_loop!(τ1, 1, ordering, α, β, τ0, ∇τ0, tags, front, κ2, dx, cs, xn, I1, Iend, LI)
    (τ1.*τ0, ordering)
end
