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
                (x0, κ², t_exact) = test_fun(hlist[i]) 
                t_pred = FEFMM.fefmm(κ², repeat([hlist[i]], ndims(κ²)), x0)
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
            κ²1D = ones(101)
            dx1D = [0.1]
            xs1D = CartesianIndex(1)
            τa1D = 0:0.1:10
            @test all(isapprox.(τa1D, FEFMM.fefmm(κ²1D, dx1D, xs1D)))
            κ²2D = ones(101,101)
            dx2D = [0.1,0.1]
            xs2D = CartesianIndex(1,1)
            τaxis = 0:0.1:10
            τa2D = sqrt.(τaxis'.^2 .+ τaxis.^2)
            @test all(isapprox.(τa2D, FEFMM.fefmm(κ²2D, dx2D, xs2D)))
            @test all(isapprox.(Float32.(τa2D), FEFMM.fefmm(Float32.(κ²2D), Float32.(dx2D), xs2D)))
            κ²3D = ones(101,101,101)
            dx3D = [0.1,0.1,0.1]
            xs3D = CartesianIndex(1,1,1)
            τa3D = zeros(101,101,101)
            for k in 1:101
                for j in 1:101
                    for i in 1:101
                        τa3D[i,j,k] = sqrt(τaxis[i]^2+τaxis[j]^2+τaxis[k]^2)
                    end
                end
            end
            @test all(isapprox.(τa3D, FEFMM.fefmm(κ²3D, dx3D, xs3D)))
            @test_broken all(isapprox.(Float32.(τa3D), FEFMM.fefmm(Float32.(κ²3D), Float32.(dx3D), xs3D)))
        end

        @testset "Convergence Tests" begin
            convergence_test([1/20, 1/40, 1/80, 1/160, 1/320], const_κ²_2D, "Constant Gradient of Squared Slowness 2D")
            convergence_test([1/20, 1/40, 1/80, 1/160, 1/320], const_v2_2D, "Constant Gradient of Squared Velocity 2D")
            convergence_test([1/20, 1/40, 1/80, 1/160, 1/320], gaussian_factor_2D, "Gaussian Factor 2D")
            convergence_test([1/20, 1/40, 1/80], const_κ²_3D, "Constant Gradient of Squared Slowness 3D")
            convergence_test([1/20, 1/40, 1/80], const_v2_3D, "Constant Gradient of Squared Velocity 3D")
            convergence_test([1/20, 1/40, 1/80], gaussian_factor_3D, "Gaussian Factor 3D")
        end
    end
end

