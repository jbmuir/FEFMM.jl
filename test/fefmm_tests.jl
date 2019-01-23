module FEFMM_Tests
    using Test
    using FEFMM

    function const_k2(h)
        #Define Medium
        a = -0.4
        s0 = 2.0
        #Set up grid
        X = 0.0:h:8.0
        Y = 0.0:h:4.0
        sx = 4.0
        sy = 0.0 
        sxi = findfirst(x->x==sx, X)
        syi = findfirst(y->y==sy, Y)
        x0 = CartesianIndex(syi, sxi)
        dist2 = ((X.-sx)').^2 .+ ((Y.-sy)).^2
        depth = 0.0*(X.-sx)'.+(Y.-sy)
        k2 = s0^2 .+ 2*a*depth
        Sb2 = s0^2 .+ a*depth
        sig2 = 2*dist2./(Sb2.+sqrt.(Sb2.^2 .- a*a*dist2))
        t_exact = Sb2.*sqrt.(sig2) .- a*a*sqrt.(sig2).^3 ./ 6
        (x0, k2, t_exact)
    end

    function const_v2(h)
        #Define Medium
        a = 1.0
        s0 = 2.0
        #Set up grid
        X = 0.0:h:8.0
        Y = 0.0:h:4.0
        sx = 4.0
        sy = 0.0 
        sxi = findfirst(x->x==sx, X)
        syi = findfirst(y->y==sy, Y)
        x0 = CartesianIndex(syi, sxi)
        dist2 = ((X.-sx)').^2 .+ ((Y.-sy)).^2
        depth = 0.0*(X.-sx)'.+(Y.-sy)
        k = 1 ./ (1/s0 .+ a*depth)

        t_exact = @. 1/a*acosh(1+1/2*s0*a*a*k*dist2)
        (x0, k.^2, t_exact)
    end

    @testset "FEFMM Tests" begin
        @testset "Solve Piecewise Quadratic" begin
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,1.0,1.0), 2.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,1.0,1.0,1.0,1.0), 1.0+1/sqrt(2))
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,1.0,100.0,1.0,1.0), 2.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(2.0,1.0,1.0,100.0,1.0), 1.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,1.0,1.0,1.0,1.0,1.0,1.0), 1.0+1/sqrt(3))
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,2.0,3.0,100.0,10.0,1.0,1.0), 2.0/3.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,2.0,3.0,10.0,100.0,1.0,1.0), 2.0/3.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,2.0,3.0,10.0,1.0,100.0,1.0), 1.0)
            @test isapprox(FEFMM.solve_piecewise_quadratic(1.0,2.0,3.0,1.0,10.0,100.0,1.0), 2.0)
        end

        @testset "Basic FEFMM Checks" begin
        κ2 = ones(101)
        dx = [0.1]
        xs = CartesianIndex(1)
        τa = 0:0.1:10
        @test all(isapprox.(τa, FEFMM.fefmm(κ2, dx, xs)[1]))
        end
    end
    #     @testset "Constant Gradient of Squared Slowness" begin
    #         h = 1/40
    #         (x0, k2, t_exact) = const_k2(h)
    #         (t_pred, ordering) = FEFMM.fefmm(k2, [h,h], x0)
    #         Profile.clear()  # in case we have any previous profiling data
    #         h = 1/160
    #         (x0, k2, t_exact) = const_k2(h)
    #         @profile FEFMM.fefmm(k2, [h,h], x0)
    #         sqrt(sum((t_pred.-t_exact).^2)/length(t_pred))
    #     end

    #     @testset "Constant Gradient of Squared Velocity" begin
    #         h = 1/40
    #         (x0, k2, t_exact) = const_v2(h)
    #         (t_pred, ordering) = fefmm(k2, [h,h], x0)
    #         sqrt(sum((t_pred.-t_exact).^2)/length(t_pred))
    #     end
    # end
end