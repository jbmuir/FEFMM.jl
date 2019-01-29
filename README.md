FEFMM.jl
=================

[![Build Status](https://travis-ci.com/jbmuir/FEFMM.jl.svg?branch=master)](https://travis-ci.com/jbmuir/FEFMM.jl)

This is a lightweight, Julia 1.0+ compatible version of Treister & Haber's Factored Eikonal Fast Marching Method. 

See: 
Eran Treister and Eldad Haber, *A fast marching algorithm for the factored eikonal equation*, Journal of Computational Physics, 324, 210-225, 2016. 

with implementation at https://github.com/JuliaInv/FactoredEikonalFastMarching.jl

This module exports a single function, e.g.

```
    using FEFMM
    k2 = ones(100,100) #slowness squared
    dx = [1.0, 1.0] # grid spacing
    x0 = CartesianIndex(1,1) # source location index
    (t, ordering) = fefmm(k2,dx,x0) #the function gives both the time estimate (t) and order in which nodes were finalized (ordering)
```
# Requirements
 - Julia 1.0+
 - DataStructures.jl
 - UnicodePlots.jl [test]
 - Test.jl [test]
 - LinearAlgebra.jl [test]
