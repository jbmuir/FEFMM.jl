module FEFMM_Tests
    using Test
    using FEFMM
    using UnicodePlots
    using LinearAlgebra: norm

    include("testmediums.jl")

    function convergence_test(hlist, test_fun, test_name)
        @testset "$test_name" begin
            linorm = zeros(length(hlist))
            l2norm = zeros(length(hlist))
            for i = 1:length(hlist)
                (x0, k2, t_exact) = test_fun(hlist[i]) 
                t_pred = FEFMM.fefmm(k2, repeat([hlist[i]], ndims(k2)), x0)
                err = t_exact.-t_pred
                linorm[i] = norm(err, Inf)
                l2norm[i] = norm(err, 2)/sqrt(length(err))
            end

            plt = lineplot(log10.(hlist), 
                        log10.(linorm),
                        title = test_name,
                        name = "L_inf norm",
                        xlabel = "log_10(h)",
                        ylabel = "log_10(Error)", 
                        xlim = [floor(minimum(log10.(hlist))), ceil(maximum(log10.(hlist)))],
                        ylim = [floor(minimum(log10.([linorm; l2norm]))), ceil(maximum(log10.([linorm; l2norm])))],
                        color = :blue)
            
            lineplot!(plt, log10.(hlist), log10.(l2norm), name = "mean L2 norm", color = :red)
            println(plt)
            @test reverse(linorm) == sort(linorm) #test for decreasing error with decreasing h
            @test reverse(l2norm) == sort(l2norm)
        end
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
            @test all(isapprox.(τa1d, FEFMM.fefmm(κ21d, dx1d, xs1d)))
            κ22d = ones(101,101)
            dx2d = [0.1,0.1]
            xs2d = CartesianIndex(1,1)
            τaxis = 0:0.1:10
            τa2d = sqrt.(τaxis'.^2 .+ τaxis.^2)
            @test all(isapprox.(τa2d, FEFMM.fefmm(κ22d, dx2d, xs2d)))
            @test all(isapprox.(Float32.(τa2d), Float32.(FEFMM.fefmm(κ22d, Float32.(dx2d), xs2d))))
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
            @test all(isapprox.(τa3d, FEFMM.fefmm(κ23d, dx3d, xs3d)))
            @test all(isapprox.(Float32.(τa3d), Float32.(FEFMM.fefmm(κ23d, Float32.(dx3d), xs3d))))
        end

        @testset "Convergence Tests" begin
            convergence_test([1/10, 1/20, 1/40, 1/80, 1/160, 1/320], const_k2_2D, "Constant Gradient of Squared Slowness 2D")
            convergence_test([1/10, 1/20, 1/40, 1/80, 1/160, 1/320], const_v2_2D, "Constant Gradient of Squared Velocity 2D")
            convergence_test([1/10, 1/20, 1/40, 1/80, 1/160, 1/320], gaussian_factor_2D, "Gaussian Factor 2D")
            convergence_test([1/10, 1/20, 1/40, 1/80], const_k2_3D, "Constant Gradient of Squared Slowness 3D")
            convergence_test([1/10, 1/20, 1/40, 1/80], const_v2_3D, "Constant Gradient of Squared Velocity 3D")
            convergence_test([1/10, 1/20, 1/40, 1/80], gaussian_factor_3D, "Gaussian Factor 3D")
        end

        @testset "N-dimensional mul_analytic! and grad_analytic" begin
            @testset "1D case" begin
                τ_1d = ones(5)
                dx1, ie1, xs1 = 0.1, 5, 1
                FEFMM.mul_analytic!(τ_1d, (dx1,), (ie1,), (xs1,))
                @test isapprox(τ_1d, [1.0, 1.1, 1.2, 1.3, 1.4])

                τ0_1d = [1.0, 1.1, 1.2, 1.3, 1.4]
                grad_1d = FEFMM.grad_analytic(τ0_1d, (dx1,), (ie1,), (xs1,))
                @test length(grad_1d) == 1
                @test size(grad_1d[1]) == (5,)
                @test isapprox(grad_1d[1], [1.0, 0.09090909, 0.16666667, 0.23076923, 0.28571429], rtol=1e-6)
            end

            @testset "2D case" begin
                τ_2d = ones(3, 3)
                dx1, dx2, ie1, ie2, xs1, xs2 = 0.1, 0.1, 3, 3, 1, 1
                FEFMM.mul_analytic!(τ_2d, (dx1, dx2), (ie1, ie2), (xs1, xs2))
                @test isapprox(τ_2d, [1.0 1.1 1.2; 1.1 1.414214 1.5; 1.2 1.5 1.732051], rtol=1e-6)

                τ0_2d = [1.0 1.1 1.2; 1.1 1.414214 1.5; 1.2 1.5 1.732051]
                grad_2d = FEFMM.grad_analytic(τ0_2d, (dx1, dx2), (ie1, ie2), (xs1, xs2))
                @test length(grad_2d) == 2
                @test all(size.(grad_2d) .== Ref((3, 3)))
                @test isapprox(grad_2d[1][1,1], 1/sqrt(2), rtol=1e-6)
                @test isapprox(grad_2d[2][1,1], 1/sqrt(2), rtol=1e-6)
            end

            @testset "3D case" begin
                τ_3d = ones(2, 2, 2)
                dx1, dx2, dx3, ie1, ie2, ie3, xs1, xs2, xs3 = 0.1, 0.1, 0.1, 2, 2, 2, 1, 1, 1
                FEFMM.mul_analytic!(τ_3d, (dx1, dx2, dx3), (ie1, ie2, ie3), (xs1, xs2, xs3))
                @test isapprox(τ_3d, cat([1.0 1.1; 1.1 1.414214], [1.1 1.414214; 1.414214 1.732051], dims=3), rtol=1e-6)

                τ0_3d = cat([1.0 1.1; 1.1 1.414214], [1.1 1.414214; 1.414214 1.732051], dims=3)
                grad_3d = FEFMM.grad_analytic(τ0_3d, (dx1, dx2, dx3), (ie1, ie2, ie3), (xs1, xs2, xs3))
                @test length(grad_3d) == 3
                @test all(size.(grad_3d) .== Ref((2, 2, 2)))
                @test isapprox(grad_3d[1][1,1,1], 1/sqrt(3), rtol=1e-6)
                @test isapprox(grad_3d[2][1,1,1], 1/sqrt(3), rtol=1e-6)
                @test isapprox(grad_3d[3][1,1,1], 1/sqrt(3), rtol=1e-6)
            end

            @testset "4D case" begin
                τ_4d = ones(2, 2, 2, 2)
                dx = (0.1, 0.1, 0.1, 0.1)
                ie = (2, 2, 2, 2)
                xs = (1, 1, 1, 1)
                FEFMM.mul_analytic!(τ_4d, dx, ie, xs)
                @test isapprox(τ_4d[2,2,2,2], sqrt(4 * 0.1^2), rtol=1e-6)

                τ0_4d = copy(τ_4d)
                grad_4d = FEFMM.grad_analytic(τ0_4d, dx, ie, xs)
                @test length(grad_4d) == 4
                @test all(size.(grad_4d) .== Ref((2, 2, 2, 2)))
                @test all(isapprox.(getindex.(grad_4d, 1, 1, 1, 1), 1/sqrt(4), rtol=1e-6))
            end
        end
    end
end

