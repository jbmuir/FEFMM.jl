function solve_node(τ₁::Array{R, N}, 
                    α::Vector{R},
                    β::Vector{R},
                    τ₀::Array{R, N},
                    ∇τ₀::Array{T},
                    tags::Array{UInt8, N},
                    κ²::Array{R, N},
                    x::S,
                    dx::Vector{<:AbstractFloat},
                    cs::Vector{S},
                    I1::S, 
                    Iend::S) where {N, R <: AbstractFloat, T <: Array{R, N}, S <: CartesianIndex{N}}
    #loop over dimensions to set \alpha and \beta ...
    for (i, s) in enumerate(cs)
        fwdok = check1fwd(τ₁, τ₀, tags, x, s, Iend)
        bwdok = check1bwd(τ₁, τ₀, tags, x, s, I1)
        if fwdok && bwdok 
            if @inbounds τ₁[x+s]*τ₀[x+s] < τ₁[x-s]*τ₀[x-s]
                bwkok = false
            else 
                fwdok = false
            end
        end
        # Note Here: working out the maths gives an additional 1/\alpha for the \beta terms
        # to reduce total operations we redefine \beta = \beta' / \alpha where \beta' is the \beta in the paper
        if fwdok 
            if check2fwd(τ₁, τ₀, tags, x, s, Iend)
                @inbounds α[i] = 1.5*τ₀[x]/dx[i]-∇τ₀[i][x]
                @inbounds β[i] = (2*τ₁[x+s]-0.5*τ₁[x+2*s])*τ₀[x]/dx[i]
            else
                @inbounds α[i] = τ₀[x]/dx[i]-∇τ₀[i][x]
                @inbounds β[i] = τ₀[x]*τ₁[x+s]/dx[i]
            end
        elseif bwdok
            if check2bwd(τ₁, τ₀, tags, x, s, I1)
                @inbounds α[i] = 1.5*τ₀[x]/dx[i]+∇τ₀[i][x]
                @inbounds β[i] = (2*τ₁[x-s]-0.5*τ₁[x-2*s])*τ₀[x]/dx[i]
            else
                @inbounds α[i] = τ₀[x]/dx[i]+∇τ₀[i][x]
                @inbounds β[i] = τ₀[x]*τ₁[x-s]/dx[i]
            end
        else
            @inbounds α[i] = zero(R)
            @inbounds β[i] = zero(R)
        end
    end
    #All of the \alpha and \beta should be set, so now we will proceed with solving the solve_quadratic
    #@inbounds solve_piecewise_quadratic(α..., β..., κ²[x])
    if N == 1
        @inbounds τ₁y = solve_piecewise_quadratic(α[1], β[1], κ²[x])
    elseif N == 2
        @inbounds τ₁y = solve_piecewise_quadratic(α[1], α[2], β[1], β[2], κ²[x])
    elseif N == 3
        @inbounds τ₁y = solve_piecewise_quadratic(α[1], α[2], α[3], β[1], β[2], β[3], κ²[x])        
    end
    τ₁y
end


function fefmm_loop!(τ₁::Array{R, N}, 
                     tags::Array{UInt8, N}, 
                     front::BinaryMinHeap{Node{R, N}},
                     α::Vector{R},
                     β::Vector{R},
                     τ₀::Array{R, N}, 
                     ∇τ₀::Array{T},
                     κ²::Array{R, N},
                     dx::Vector{R},
                     cs::Vector{S},
                     xn::Vector{S},
                     I1::S, 
                     Iend::S) where {N, R <: AbstractFloat, T <: Array{R, N}, S <: CartesianIndex{N}}
    while isempty(front) == false
        x = pop!(front).ind
        @inbounds if tags[x] < 0x3# as we are using a non-mutable minheap, we might have already set this node to known on another pass (i.e. this could be an old worthless version of the node). This check makes sure we don't accidentally override something we already have computed
            @inbounds tags[x] = 0x3
            set_neighbours!(xn, x, tags, cs, I1, Iend)
            for y in xn 
                if y != x # to avoid allocating neighbours repeatedly we reuse xn but have y = x as the "neighbour not acceptable" flag value
                    @inbounds tags[y] = 0x2 #move neighbours to front
                    τ₁y = solve_node(τ₁, α, β, τ₀, ∇τ₀, tags, κ², y, dx, cs, I1, Iend)
                    @inbounds if τ₁y < τ₁[y]
                            #only update and add new node to heap if the new solution is smaller 
                        @inbounds τ₁[y] = τ₁y
                        @inbounds push!(front, Node(y, τ₁[y]*τ₀[y]))
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
τ₁τ₀ = travel time array from source

Tags: 
1 = unknown
2 = front
3 = known
"""

function fefmm(κ²::Array{R, N},
               dx::Vector{<:AbstractFloat},
               xs::CartesianIndex{N}) where {R <: AbstractFloat, N}
    #initialization
    #main working arrays
    τ₀ = ones(R, size(κ²)) # analytic solution array
    τ₁ = R(Inf).*ones(R, size(κ²)) # correction array
    tags = ones(UInt8, size(κ²)) # node bookkeeping array

    #derive indices and preallocate small working arrays
    cs = cartstrides(κ²)
    inds = CartesianIndices(κ²)
    I1 = first(inds)
    Iend = last(inds)
    α = Array{R}(undef, ndims(τ₁))
    β = similar(α)
    xn = Array{typeof(I1)}(undef, N*2)

    #setup main working arrays to starting configuration
    mul_analytic!(τ₀, dx..., size(κ²)..., Tuple(xs)...)
    ∇τ₀ = grad_analytic(τ₀, dx..., size(κ²)..., Tuple(xs)...)
    τ₁[xs] = sqrt(κ²[xs])
    tags[xs] = 0x2 
    front = BinaryMinHeap{Node{R,N}}()
    push!(front, Node(inds[xs], τ₁[xs]*τ₀[xs])) 

    #main loop
    fefmm_loop!(τ₁, tags, front, α, β, τ₀, ∇τ₀, κ², dx, cs, xn, I1, Iend)
    τ₁.*τ₀
end

