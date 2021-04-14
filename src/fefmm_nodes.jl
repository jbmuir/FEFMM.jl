struct Node{R <: AbstractFloat, N}
    ind::CartesianIndex{N}
    val::R
end

isless(a::Node, b::Node) = isless(a.val, b.val)

