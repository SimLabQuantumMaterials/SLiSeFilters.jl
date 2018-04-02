"F, the residual level Fct"
function residual( z::Array{T,1}, alpha::Array{T,1}, c::SLiSeFilters.Config ):T
    X = calcX(z, c);
    W = calcW(z, c);
    Z = calcZ(z, c);
    Y = calcY(z, c);
    theta = calctheta(z, c);

    XZa = (X - Z) * alpha
    WYa = (W - Y) * conj(alpha)

    return (2*real(alpha' * (XZa + WYa -2*theta)) + real(caseh_h( c )))[1]
end

"The gradient of F"
function Fgradient( z::Array{T,1}, alpha::Array{T,1}, c::SLiSeFilters.Config )#:Tuple{Array{T,1},Array{T,1}}
    X = calcX(z, c);
    Z = calcZ(z, c);
    theta = calctheta(z, c);
    Wc = calcW(conj(z), c);
    Yc = calcY(conj(z), c);

    da = 4*(alpha' * (X - Z) + alpha.' * (Wc - Yc) - theta' );

    dYc = calcDeltaYc( z, c );
    dWc = calcDeltaWc( z, c );
    dX = calcDeltaX( z, c );
    dZ = calcDeltaZ( z, c );
    dTheta = calcDeltaTheta( z, c );

    dz = 4*(alpha' * (dX - dZ) + alpha.' * (dWc - dYc) - dTheta') .* alpha.'; # elementwise multiplication

    # we conjugate as we want column vectors
    return Array{T,1}(dz'[:,1]),Array{T,1}(da'[:,1])
end
