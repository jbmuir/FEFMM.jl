FEFMM.jl
=================

This is a lightweight, Julia 1.0+ compatible version of Treister & Haber's Factored Eikonal Fast Marching Method. 

See Treister & Haber "A Fast marching algorithm for the factored eikonal equation" https://arxiv.org/abs/1607.00973 with implementation at https://github.com/JuliaInv/FactoredEikonalFastMarching.jl

This module exports a single method, e.g.

```
    using FEFMM
    k2 = ones(100,100) #slowness squared
    dx = [1.0, 1.0] # grid spacing
    x0 = CartesianIndex(1,1) # source location index
    t = fefmm(k2,dx,x0)
```
# Requirements
 - Julia 0.7+
 - DataStructures.jl
 - Plots.jl [test]
 - Test.jl [test]
 - LinearAlgebra.jl [test]