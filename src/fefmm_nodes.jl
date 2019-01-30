struct Node{R <: AbstractFloat, N}
    ind::CartesianIndex{N}
    val::R
end

<(a::Node, b::Node) = a.val < b.val
