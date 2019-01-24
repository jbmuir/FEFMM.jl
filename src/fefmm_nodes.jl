struct Node{R <: Real, N}
    ind::CartesianIndex{N}
    val::R
end

<(a::Node, b::Node) = a.val < b.val
