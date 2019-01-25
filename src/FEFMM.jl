module FEFMM

    import Base: < 
    using DataStructures
    
    export fefmm

    #small helper functions
    include("fefmm_nodes.jl")
    include("fefmm_neighbour_checks.jl")
    include("fefmm_piecewise_quadratic.jl")
    include("fefmm_utilities.jl")
    #core functions (initialization, main loop, node solve)
    include("factoredeikonalfmm.jl")

end 
