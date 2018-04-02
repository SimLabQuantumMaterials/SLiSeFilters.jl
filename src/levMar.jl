function levMar(z::Array{T,1}, alpha::Array{T,1}, c)
    n :: UInt64 = length(z);
    const maxiter :: UInt64 = 1000000;
    mu :: Float64 = 200;
    v :: UInt64 = 2;

    H = Array{T}(8*n,8*n);

    resid_old = residual(z, alpha, c);
    H,d =  hessianApprox(z, alpha, c );

    for j = 1:maxiter

        alpha_old = alpha
        z_old = z

        if isinf( mu ) || isnan( mu )
            @printf "mu too large or too small: %e\n" mu
            break;
        end

        s = (im * imag(H) - diagm(diag(H))*mu)\d;

        z_cand = z + s[1:length(z)];
        alpha_cand = alpha + s[length(z)+1:2*length(z)];
        resid_cand = residual(z_cand, alpha_cand, c);

        # this should be positive
        #LL =  real(dot(s , ( mu * s + d )));

        gain = (resid_old - resid_cand)

        if gain[1] > 0
            v = 2;

            mu = mu * 1/2;
            z = z_cand;
            alpha = alpha_cand;
            resid_old = residual(z, alpha, c);

            H,d =  hessianApprox( z, alpha, c );

        else
            mu = mu * v; v = 2*v;
        end

        if mod(j,10000) == 0
            @printf "iter: %d\tresid: %1.15e\tmu: %e\tv: %e\n" j resid_old[1] mu v
        end
    end

    @printf "resid: %e\tmu: %e\tv: %e\n" resid_old[1] mu v
    return z, alpha;
end


# the approximation to the hessian required for levMar
function hessianApprox(z, a, c)

    ddX = calcddX( z, c );
    ddZ = calcddZ( z, c );

    X = calcX( z, c);
    Z = calcZ( z, c);

    dX = calcDeltaX( z, c );
    dZ = calcDeltaZ( z, c );

    H2 = 2*[  (diagm(a)' * ( ddX - ddZ) * diagm(a)) (diagm(a)' * (dX - dZ)')  ;
              ( ( dX -  dZ ) * diagm(a) )          (X - Z)                    ];

    dz, da = Fgradient( z, a, c );

    H = H2;
    d = (0.5.*[dz;da]);
    d = (d[:])

    return H, d;
end
