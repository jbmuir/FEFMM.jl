module FEFMM

    import Base: < 
    using DataStructures: isempty, binary_minheap, BinaryHeap, LessThan
    
    export fefmm

    #small helper functions
    include("fefmmutilities.jl")
    #core functions (initialization, main loop, node solve)
    include("factoredeikonalfmm.jl")

end 
