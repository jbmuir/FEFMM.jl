function solve_node(τ1::Array{R}, 
                          α::Array{R},
                          β::Array{R},
                          τ0::Array{R},
                          ∇τ0::Array{T},
                          tags::Array{UInt8},
                          κ2::Array{R},
                          x::S,
                          dx::Array{R}, 
                          cs::Array{S}, 
                          I1::S, 
                          Iend::S) where {R <: Real, T <: Array{R}, S <: CartesianIndex}
    #loop over dimensions to set \alpha and \beta ...
    l = length(tags)
    for (i, s) in enumerate(cs)
        fwdok, fwdval = check1fwd(τ1, τ0, tags, x, s, Iend)
        bwdok, bwdval = check1bwd(τ1, τ0, tags, x, s, I1)
        if fwdok && bwdok 
            if fwdval < bwdval
                bwkok = false
            else 
                fwdok = false
            end
        end
        if fwdok 
            if check2fwd(τ1, τ0, tags, x, s, Iend)
                α[i] = 1.5*τ0[x]/dx[i]-∇τ0[i][x]
                β[i] = (2*τ1[x+s]-0.5*τ1[x+2*s])*τ0[x]/(dx[i]*α[i])
            else
                α[i] = τ0[x]/dx[i]-∇τ0[i][x]
                β[i] = τ0[x]*τ1[x+s]/(dx[i]*α[i])
            end
        elseif bwdok
            if check2bwd(τ1, τ0, tags, x, s, I1)
                α[i] = 1.5*τ0[x]/dx[i]+∇τ0[i][x]
                β[i] = (2*τ1[x-s]-0.5*τ1[x-2*s])*τ0[x]/(dx[i]*α[i])
            else
                α[i] = τ0[x]/dx[i]+∇τ0[i][x]
                β[i] = τ0[x]*τ1[x-s]/(dx[i]*α[i])
            end
        else
            α[i] = zero(R)
            β[i] = zero(R)
        end
    end
    #All of the \alpha and \beta should be set, so now we will proceed with solving the solve_quadratic
    solve_piecewise_quadratic(α, β, κ2[x])
end


function fefmm_loop!(τ1::Array{R}, 
                     ordering::Array{Int},
                     τ0::Array{R}, 
                     ∇τ0::Array{T},
                     tags::Array{UInt8}, 
                     front::BinaryHeap{Node{R}, LessThan}, 
                     κ2::Array{R},
                     dx::Array{R}, 
                     cs::Array{S}, 
                     I1::S, 
                     Iend::S) where {R <: Real, T <: Array{R}, S <: CartesianIndex}
    α = Array{R}(undef,length(size(τ1)))
    β = similar(α)
    LI = LinearIndices(τ1)
    while isempty(front) == false
        x = pop!(front).ind
        if tags[x] < 0x3 # as we are using a non-mutable minheap, we might have already set this node to known on another pass (i.e. this could be an old worthless version of the node). This check makes sure we don't accidentally override something we already have computed
            push!(ordering, LI[x])
            tags[x] = 0x3
            xn = neighbours(x, tags, cs, I1, Iend)
            for y in xn 
                tags[y] = 0x2 #move neighbours to front
                τ1y = solve_node(τ1, α, β, τ0, ∇τ0, tags, κ2, y, dx, cs, I1, Iend)
                if τ1y < τ1[y]
                        #only update and add new node to heap if the new solution is smaller 
                    τ1[y] = τ1y
                    push!(front, Node(y, τ1[y]*τ0[y]))
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

function fefmm(κ2::Array{R}, dx::Array{R}, x0::CartesianIndex) where {R <: Real}
    #initialization
    ordering = Array{Int}(undef,0)
    cs = cartstrides(κ2)
    inds = CartesianIndices(κ2)
    I1 = first(inds)
    Iend = last(inds)
    τ0 = ones(R, size(κ2))
    mul_analytic!(τ0, dx, x0)
    ∇τ0 = Array{Array{Float64,length(size(κ2))}}(undef, 0)
    for i in 1:length(size(τ0))
        ∇τ0i = ones(R, size(κ2))
        mul_grad_analytic!(∇τ0i, dx, x0, i)
        push!(∇τ0, ∇τ0i)
    end
    τ1 = Inf.*ones(R, size(κ2))
    τ1[x0] = κ2[x0]
    tags = ones(UInt8,size(κ2))
    tags[x0] = 0x2
    front = binary_minheap(Node{R})
    push!(front, Node(inds[x0], τ1[x0]*τ0[x0])) 
    #main loop
    fefmm_loop!(τ1, ordering, τ0, ∇τ0, tags, front, κ2, dx, cs, I1, Iend)
    τ1
    ordermat = ones(size(κ2))
    for (i, x) in enumerate(ordering)
        ordermat[x] = i
    end
    (τ1.*τ0, ordering)
end

κ2 = ones(5,5)
R = Float64
dx = [1.0, 1.0]
x0 = CartesianIndex(1,1)