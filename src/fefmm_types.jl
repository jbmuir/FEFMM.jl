struct Node{R <: AbstractFloat, N}
    ind::CartesianIndex{N}
    val::R
end

isless(a::Node, b::Node) = isless(a.val, b.val)

abstract type SlownessGrid end

struct CartesianGrid1D{R <: AbstractFloat} <: SlownessGrid
    
end

struct CartesianGrid2D{R <: AbstractFloat} <: SlownessGrid
end

struct CartesianGrid3D{R <: AbstractFloat} <: SlownessGrid
end

struct SphericalGrid2D{R <: AbstractFloat} <: SlownessGrid
end

struct SphericalGrid3D{R <: AbstractFloat} <: SlownessGrid
end
