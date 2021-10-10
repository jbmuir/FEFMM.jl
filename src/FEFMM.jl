module FEFMM

    import Base: isless
    using DataStructures
    
    export fefmm
    export CartesianGrid1D, 
           CartesianGrid2D, 
           CartesianGrid3D, 
           SphericalGrid2D, 
           SphericalGrid3D

    #small helper functions
    include("fefmm_types.jl")
    include("fefmm_neighbour_checks.jl")
    include("fefmm_piecewise_quadratic.jl")
    include("fefmm_utilities.jl")
    #core functions (initialization, main loop, node solve)
    include("factoredeikonalfmm.jl")

end 
