# Matrices and Vectors for the Gradients

function calctheta( z::Array{T,1}, c::SLiSeFilters.Config )
    t = zeros(Complex{Float64},length(z),1);
    for k=1:length(z)
        t[k] = theta_k( z[k], c );
    end
    return t
end

function calcthetaP( z, c::SLiSeFilters.Config )
    t = zeros(Complex{Float64},length(z),1);
    for k=1:length(z)
        t[k] = thetaP_k( z[k], c::SLiSeFilters.Config );
    end
    return t
end


# z is Hermitian
function calcZ( z, c )
    Z = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            Z[k,l] = Z_kl( z[k], z[l], c );
        end
    end
    Z = Z + tril(Z, -1)';
end

# y is real symmetric
function calcY( z, c::SLiSeFilters.Config )
    Y = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            Y[k,l] = Y_kl( z[k], z[l], c );
        end
    end
    Y = Y + transpose(tril(Y, -1));
end

# W is real symmetric
function  calcW( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            if k == l
                W[k,l] = W_kk( z[k], c );
            else
                W[k,l] = W_kl( z[k], z[l], c );
            end
        end
    end
    W = W + transpose(tril(W, -1));
end

# W is real symmetric
function calcWtp( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            if k == l
                W[k,l] = Wtp_kk( z[k], c );
            else
                W[k,l] = Wtp_kl( z[k], z[l], c );
            end
        end
    end
    W = W + transpose(tril(W, -1));
end

# X is hermitian
function calcX( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            X[k,l] = X_kl( z[k], z[l], c );
        end
    end
    X = X + tril(X, -1)';
end

# Xtp is hermitian
function calcXtp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z), length(z));
    for k=1:length(z)
        for l=1:k
            X[k,l] = Xtp_kl( z[k], z[l], c );
        end
    end
    X = X + tril(X, -1)';
end


function theta_k( z1, c::SLiSeFilters.Config )
    case1( conj(z1), c );
end

function thetaP_k( z1, c::SLiSeFilters.Config )
    x = case1( -conj(z1), c );
end

function Wtp_kk( z1, c::SLiSeFilters.Config )
    case2( -z1, c );
end

function W_kk( z1, c::SLiSeFilters.Config )
    case2( conj(z1), c );
end

function Z_kl( z1, z2, c::SLiSeFilters.Config )
    case3( conj(z1), -z2, c );
end

function  Y_kl( z1, z2, c::SLiSeFilters.Config )
    case3( conj(z1), -conj(z2), c );
end

# 1 / (t-conj(z1)) / (t-conj(z2)) -> 1 / (t-z1) / (t-z2)
function  W_kl( z1, z2, c::SLiSeFilters.Config )
    x = case3( conj(z1), conj(z2), c );
end

# 1 / (t-conj(z1)) / (t-conj(z2)) -> 1 / (t-z1) / (t-z2)
function  Wtp_kl( z1, z2, c::SLiSeFilters.Config )
    x = case3( -z1, -z2, c );
end

#        % 1 / (t-conj(z1)) / (t-z2) -> 1 / (t-z1) / (t-z2)
function  X_kl( z1, z2, c::SLiSeFilters.Config )
    x = case3( conj(z1), z2, c );
end

#        % 1 / (t+z1)) / (t+conj(z2)) -> 1 / (t-z1) / (t-z2)
function Xtp_kl( z1, z2, c::SLiSeFilters.Config )
    x = case3( -z1, -conj(z2), c );
end


############################################################
function calcddW( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z),length(z));
    for k = 1:length(z)
        for l = 1:length(z)
            if z[k] == z[l]
                W[k,k] = case7( z[k], c );
            else
                W[k,l] = case6( z[k], z[l], c );
            end
        end
    end
    return W
end

function calcddX( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z),length(z));
    for k = 1:length(z)
        for l = 1:length(z)
            zk = conj(z[k]);
            zl = z[l];
            if zk == zl
                W[k,l] = case7( zk, c );
            else
                W[k,l] = case6( zk, zl, c );
            end
        end
    end
    return W
end

function calcddY( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z),length(z));
    for k = 1:length(z)
        for l = 1:length(z)
            zk = z[k];
            zl = -z[l];
            if zk == zl
                W[k,l] = -case7( zk, c );
            else
                W[k,l] = -case6( zk, zl, c );
            end
        end
    end
    return W
end

function calcddZ( z, c::SLiSeFilters.Config )
    W = zeros(Complex{Float64},length(z),length(z));
    for k = 1:length(z)
        for l = 1:length(z)
            zk = conj(z[k]);
            zl = -z[l];
            if zk == zl
                W[k,l] = -case7( zk, c );
            else
                W[k,l] = -case6( zk, zl, c );
            end
        end
    end
    return W
end

function calcDeltaX( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = case4( z[l], conj(z[k]), c );
        end
    end
    return X
end


function calcDeltaXp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = -case4( -z[l], -conj(z[k]), c );
        end
    end
    return X
end

function calcDeltaYc( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = -case4( -z[l], z[k], c );
        end
    end
    return X
end

function calcDeltaYcp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = case4( z[l], -z[k], c );
        end
    end
    return X
end


function calcDeltaWc( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            if k == l
                X[k,l] = case5( z[l], c );
            else
                X[k,l] = case4( z[l], z[k], c );
            end
        end
    end
    return X
end


function calcDeltaWcp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            if k == l
                X[k,l] = -case5( -z[l], c );
            else
                X[k,l] = -case4( -z[l], -z[k], c );
            end
        end
    end
    return X
end

function calcDeltaZ( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = -case4( -z[l], conj(z[k]), c );
        end
    end
    return X
end

function calcDeltaZp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),length(z));
    for k=1:length(z)
        for l=1:length(z)
            X[k,l] = case4( z[l], -conj(z[k]), c );
        end
    end
    return X
end

function calcDeltaThetac( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),1);
    for k=1:length(z)
        X[k] = case2_h( z[k], c );
    end
    return X
end

function calcDeltaTheta( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),1);
    for k=1:length(z)
        X[k] = case2_h( conj(z[k]), c );
    end
    return X
end

function calcDeltaThetacp( z, c::SLiSeFilters.Config )
    X = zeros(Complex{Float64},length(z),1);
    for k=1:length(z)
        X[k] = -case2_h( -z[k], c );
    end
    return X
end
