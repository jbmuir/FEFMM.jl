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
        κ21d = ones(101)
        dx1d = [0.1]
        xs1d = CartesianIndex(1)
        τa1d = 0:0.1:10
        @test all(isapprox.(τa1d, FEFMM.fefmm(κ21d, dx1d, xs1d)[1]))
        κ22d = ones(101,101)
        dx2d = [0.1,0.1]
        xs2d = CartesianIndex(1,1)
        τaxis = 0:0.1:10
        τa2d = sqrt.(τaxis'.^2 .+ τaxis.^2)
        @test all(isapprox.(τa2d, FEFMM.fefmm(κ22d, dx2d, xs2d)[1]))
        κ23d = ones(101,101,101)
        dx3d = [0.1,0.1,0.1]
        xs3d = CartesianIndex(1,1,1)
        τa3d = zeros(101,101,101)
        for k in 1:101
            for j in 1:101
                for i in 1:101
                    τa3d[i,j,k] = sqrt(τaxis[i]^2+τaxis[j]^2+τaxis[k]^2)
                end
            end
        end
        @test all(isapprox.(τa3d, FEFMM.fefmm(κ23d, dx3d, xs3d)[1]))
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