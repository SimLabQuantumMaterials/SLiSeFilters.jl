"F, the residual level Fct"
function residual( z::Array{T,1}, alpha::Array{T,1}, c::SLiSeFilters.Config ):T
    X = calcX(z, c);
    W = calcW(z, c);
    Z = calcZ(z, c);
    Y = calcY(z, c);
    theta = calctheta(z, c);
    
    XZa = (X - Z) * alpha
    WYa = (W - Y) * conj(alpha)

    return (2*real(alpha' * (XZa + WYa -2*theta)) .+ real(caseh_h( c )))[1]
end

"The gradient of F"
function Fgradient( z::Array{T,1}, alpha::Array{T,1}, c::SLiSeFilters.Config )
    X = calcX(z, c);
    Z = calcZ(z, c);
    theta = calctheta(z, c);
    Wc = calcW(conj(z), c);
    Yc = calcY(conj(z), c);
    
    da = 4*(alpha' * (X - Z) + transpose(alpha) * (Wc - Yc) - theta' );
    
    dYc = calcDeltaYc( z, c );
    dWc = calcDeltaWc( z, c );
    dX = calcDeltaX( z, c );
    dZ = calcDeltaZ( z, c );
    dTheta = calcDeltaTheta( z, c );
    
    dz = 4*(alpha' * (dX - dZ) + transpose(alpha) * (dWc - dYc) - dTheta') .* transpose(alpha); # elementwise multiplication
    
    # we conjugate as we want column vectors
    return dz',da'
end

# the approximation to the hessian required for levMar
function hessianApprox( z::Array{T,1}, a::Array{T,1}, c::SLiSeFilters.Config )
    
    ddX = calcddX( z, c );
    ddZ = calcddZ( z, c );
    
    X = calcX( z, c);
    Z = calcZ( z, c);
    
    dX = calcDeltaX( z, c );
    dZ = calcDeltaZ( z, c );
    
    H2 = 2*[  (diagm(a)' * ( ddX - ddZ) * diagm(a)) (diagm(a)' * (dX - dZ)')  ;
              ( ( dX -  dZ ) * diagm(a) )          (X - Z)                  ];

    dz, da = Fgradient( z, a, c );
    
    H = H2;
    d = (0.5.*[dz;da]);
    d = (d[:])
    
    return H, d;
end
