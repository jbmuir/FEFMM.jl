module FEFMM

    import Base: < 
    using DataStructures
    
    export fefmm

    #small helper functions
    include("fefmmutilities.jl")
    #core functions (initialization, main loop, node solve)
    include("factoredeikonalfmm.jl")

end 
