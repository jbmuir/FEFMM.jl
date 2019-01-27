function const_k2_2D(h)
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
    k2 = s0^2 .+ 2*a*depth
    Sb2 = s0^2 .+ a*depth
    sig2 = 2*dist2./(Sb2.+sqrt.(Sb2.^2 .- a*a*dist2))
    t_exact = Sb2.*sqrt.(sig2) .- a*a*sqrt.(sig2).^3 ./ 6
    (x0, k2, t_exact)
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
    t_exact = @. 1/a*acosh(1+1/2*s0*a*a*k*dist2)
    (x0, k.^2, t_exact)
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
    t1_exact = @. exp(-(sigmaz*((Z-x1z)^2)+sigmax*((X'-x1x)^2)))/2+1/2
    t1_dz_exact = @. -2*sigmaz*(Z-x1z)*(t1_exact-1/2)
    t1_dx_exact = @. -2*sigmax*(X'-x1x)*(t1_exact-1/2)

    t0 = ones(size(t1_exact))
    FEFMM.mul_analytic!(t0, h, h, size(t1_exact)..., szi, sxi)
    t0_dz, t0_dx = FEFMM.grad_analytic(t0, h, h, size(t1_exact)..., szi, sxi)

    t_dz = @. t0*t1_dz_exact+t0_dz*t1_exact
    t_dx = @. t0*t1_dx_exact+t0_dx*t1_exact
    k2 = @. t_dz*t_dz+t_dx*t_dx
    (x0, k2, t1_exact.*t0)
end

function const_k2_3D(h)
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

    k2 = @. s0^2+2*a*depth
    Sb2 = @. s0^2+a*depth
    sig2 = @. 2*dist2/(Sb2+sqrt(Sb2^2-a*a*dist2))
    t_exact = @. Sb2*sqrt(sig2)-a*a*(sqrt(sig2)^3)/6
    (x0, k2, t_exact)
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
    t_exact = @. 1/a*acosh(1+1/2*s0*a*a*k*dist2)
    (x0, k.^2, t_exact)
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

    t1_exact = zeros(length(X1), length(X2), length(X3))
    t1_dx1_exact = zeros(length(X1), length(X2), length(X3))
    t1_dx2_exact = zeros(length(X1), length(X2), length(X3))
    t1_dx3_exact = zeros(length(X1), length(X2), length(X3))

    # get medium
    for k = 1:length(X3)
        for j = 1:length(X2)
            for i = 1:length(X1)
                t1_exact[i,j,k] = exp(-(sigmax1*(X1[i]-x11)^2+sigmax2*(X2[j]-x12)^2+sigmax3*(X3[k]-x13)^2))/2+1/2
                t1_dx1_exact[i,j,k] = -2*sigmax1*(X1[i]-x11)*(t1_exact[i,j,k]-1/2)
                t1_dx2_exact[i,j,k] = -2*sigmax2*(X2[j]-x12)*(t1_exact[i,j,k]-1/2)
                t1_dx3_exact[i,j,k] = -2*sigmax3*(X3[k]-x13)*(t1_exact[i,j,k]-1/2)
            end
        end
    end

    t0 = ones(size(t1_exact))
    FEFMM.mul_analytic!(t0, h, h, h, size(t1_exact)..., sx1i, sx2i, sx3i)
    t0_dx1, t0_dx2, t0_dx3 = FEFMM.grad_analytic(t0, h, h, h, size(t1_exact)...,  sx1i, sx2i, sx3i)

    t_dx1 = @. t0*t1_dx1_exact+t0_dx1*t1_exact
    t_dx2 = @. t0*t1_dx2_exact+t0_dx2*t1_exact
    t_dx3 = @. t0*t1_dx3_exact+t0_dx3*t1_exact

    k2 = @. t_dx1*t_dx1+t_dx2*t_dx2+t_dx3*t_dx3
    (x0, k2, t1_exact.*t0)
end