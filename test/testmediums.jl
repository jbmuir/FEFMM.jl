function const_κ²_2D(h)
    #Define Medium
    a = -0.4
    s0 = 2.0
    #Set up grid
    X = 0.0:h:8.0
    Z = 0.0:h:4.0
    sx = 4.0
    sz = 0.0 
    sxi = findfirst(x->x==sx, X)
    szi = findfirst(z->z==sz, Z)
    x0 = CartesianIndex(szi, sxi)
    dist2 = ((X.-sx)').^2 .+ ((Z.-sz)).^2
    depth = 0.0*(X.-sx)'.+(Z.-sz)
    κ² = s0^2 .+ 2*a*depth
    Sb2 = s0^2 .+ a*depth
    sig2 = 2*dist2./(Sb2.+sqrt.(Sb2.^2 .- a*a*dist2))
    τ_exact = Sb2.*sqrt.(sig2) .- a*a*sqrt.(sig2).^3 ./ 6
    (x0, κ², τ_exact)
end

function const_v2_2D(h)
    #Define Medium
    a = 1.0
    s0 = 2.0
    #Set up grid
    X = 0.0:h:8.0
    Z = 0.0:h:4.0
    sx = 4.0
    sz = 0.0 
    sxi = findfirst(x->x==sx, X)
    szi = findfirst(z->z==sz, Z)
    x0 = CartesianIndex(szi, sxi)
    dist2 = ((X.-sx)').^2 .+ ((Z.-sz)).^2
    depth = 0.0*(X.-sx)'.+(Z.-sz)
    k = 1 ./ (1/s0 .+ a*depth)
    τ_exact = @. 1/a*acosh(1+1/2*s0*a*a*k*dist2)
    (x0, k.^2, τ_exact)
end

function gaussian_factor_2D(h)
    sigmax = 0.4
    sigmaz = 0.1
    x1x_attempt = 2
    x1z_attempt = 4/3
    sx = 2.0
    sz = 1.0
    X = 0.0:h:8.0
    Z = 0.0:h:4.0
    sxi = findfirst(x->x==sx, X)#-1 # -1 here is for consistency with Treister & Haber paper / implementation
    szi = findfirst(z->z==sz, Z)#-1
    x1x = X[argmin(abs.(X.-x1x_attempt))]#-h # -h here is for consistency with Treister & Haber paper / implementation
    x1z = Z[argmin(abs.(Z.-x1z_attempt))]#-h

    x0 = CartesianIndex(szi, sxi)
    τ₁_exact = @. exp(-(sigmaz*((Z-x1z)^2)+sigmax*((X'-x1x)^2)))/2+1/2
    τ₁_dz_exact = @. -2*sigmaz*(Z-x1z)*(τ₁_exact-1/2)
    τ₁_dx_exact = @. -2*sigmax*(X'-x1x)*(τ₁_exact-1/2)

    τ₀ = ones(size(τ₁_exact))
    FEFMM.mul_analytic!(τ₀, h, h, size(τ₁_exact)..., szi, sxi)
    τ₀_dz, τ₀_dx = FEFMM.grad_analytic(τ₀, h, h, size(τ₁_exact)..., szi, sxi)

    τ_dz = @. τ₀*τ₁_dz_exact+τ₀_dz*τ₁_exact
    τ_dx = @. τ₀*τ₁_dx_exact+τ₀_dx*τ₁_exact
    κ² = @. τ_dz*τ_dz+τ_dx*τ_dx
    (x0, κ², τ₁_exact.*τ₀)
end

function const_κ²_3D(h)
    #Define Medium
    a = -1.65
    s0 = 2.0
    #Set up grid
    X1 = 0.0:h:0.8
    X2 = 0.0:h:1.6
    X3 = 0.0:h:1.6
    sx1 = 0.0
    sx2 = 0.8
    sx3 = 0.8
    sx1i = findfirst(x->x==sx1, X1)
    sx2i = findfirst(x->x==sx2, X2)
    sx3i = findfirst(x->x==sx3, X3)
    x0 = CartesianIndex(sx1i, sx2i, sx3i)

    dist2 = zeros(length(X1), length(X2), length(X3))
    depth = zeros(length(X1), length(X2), length(X3))
    for k = 1:length(X3)
        for j = 1:length(X2)
            for i = 1:length(X1)
                dist2[i,j,k] = (X1[i]-sx1)^2+(X2[j]-sx2)^2+(X3[k]-sx3)^2
                depth[i,j,k] = (X1[i]-sx1)
            end
        end
    end

    κ² = @. s0^2+2*a*depth
    Sb2 = @. s0^2+a*depth
    sig2 = @. 2*dist2/(Sb2+sqrt(Sb2^2-a*a*dist2))
    τ_exact = @. Sb2*sqrt(sig2)-a*a*(sqrt(sig2)^3)/6
    (x0, κ², τ_exact)
end

function const_v2_3D(h)
    #Define Medium
    a = 1.0
    s0 = 2.0
    #Set up grid
    X1 = 0.0:h:0.8
    X2 = 0.0:h:1.6
    X3 = 0.0:h:1.6
    sx1 = 0.0
    sx2 = 0.8
    sx3 = 0.8
    sx1i = findfirst(x->x==sx1, X1)
    sx2i = findfirst(x->x==sx2, X2)
    sx3i = findfirst(x->x==sx3, X3)
    x0 = CartesianIndex(sx1i, sx2i, sx3i)

    dist2 = zeros(length(X1), length(X2), length(X3))
    depth = zeros(length(X1), length(X2), length(X3))
    for k = 1:length(X3)
        for j = 1:length(X2)
            for i = 1:length(X1)
                dist2[i,j,k] = (X1[i]-sx1)^2+(X2[j]-sx2)^2+(X3[k]-sx3)^2
                depth[i,j,k] = (X1[i]-sx1)
            end
        end
    end
    
    #get parameters
    k = @. 1/(1/s0+a*depth)
    τ_exact = @. 1/a*acosh(1+1/2*s0*a*a*k*dist2)
    (x0, k.^2, τ_exact)
end

function gaussian_factor_3D(h)
    #define medium
    sigmax1 = 0.2
    sigmax2 = 0.4
    sigmax3 = 0.1
    x11_attempt = 0.4
    x12_attempt = 1.6/3
    x13_attempt = 0.4
    sx1 = 0.2
    sx2 = 0.4
    sx3 = 0.4
    #set up grid
    X1 = 0.0:h:0.8
    X2 = 0.0:h:1.6
    X3 = 0.0:h:1.6
    sx1i = findfirst(x->x==sx1, X1)
    sx2i = findfirst(x->x==sx2, X2)
    sx3i = findfirst(x->x==sx3, X3)
    x0 = CartesianIndex(sx1i, sx2i, sx3i)

    x11 = X1[argmin(abs.(X2.-x11_attempt))]#-h # -h here is for consistency with Treister & Haber paper / implementation
    x12 = X2[argmin(abs.(X2.-x12_attempt))]#-h
    x13 = X3[argmin(abs.(X3.-x13_attempt))]

    τ₁_exact = zeros(length(X1), length(X2), length(X3))
    τ₁_dx1_exact = zeros(length(X1), length(X2), length(X3))
    τ₁_dx2_exact = zeros(length(X1), length(X2), length(X3))
    τ₁_dx3_exact = zeros(length(X1), length(X2), length(X3))

    # get medium
    for k = 1:length(X3)
        for j = 1:length(X2)
            for i = 1:length(X1)
                τ₁_exact[i,j,k] = exp(-(sigmax1*(X1[i]-x11)^2+sigmax2*(X2[j]-x12)^2+sigmax3*(X3[k]-x13)^2))/2+1/2
                τ₁_dx1_exact[i,j,k] = -2*sigmax1*(X1[i]-x11)*(τ₁_exact[i,j,k]-1/2)
                τ₁_dx2_exact[i,j,k] = -2*sigmax2*(X2[j]-x12)*(τ₁_exact[i,j,k]-1/2)
                τ₁_dx3_exact[i,j,k] = -2*sigmax3*(X3[k]-x13)*(τ₁_exact[i,j,k]-1/2)
            end
        end
    end

    τ₀ = ones(size(τ₁_exact))
    FEFMM.mul_analytic!(τ₀, h, h, h, size(τ₁_exact)..., sx1i, sx2i, sx3i)
    τ₀_dx1, τ₀_dx2, τ₀_dx3 = FEFMM.grad_analytic(τ₀, h, h, h, size(τ₁_exact)...,  sx1i, sx2i, sx3i)

    τ_dx1 = @. τ₀*τ₁_dx1_exact+τ₀_dx1*τ₁_exact
    τ_dx2 = @. τ₀*τ₁_dx2_exact+τ₀_dx2*τ₁_exact
    τ_dx3 = @. τ₀*τ₁_dx3_exact+τ₀_dx3*τ₁_exact

    κ² = @. τ_dx1*τ_dx1+τ_dx2*τ_dx2+τ_dx3*τ_dx3
    (x0, κ², τ₁_exact.*τ₀)
end

