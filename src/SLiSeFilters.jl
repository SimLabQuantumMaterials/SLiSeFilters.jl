module SLiSeFilters

using Printf, LinearAlgebra

export createConfig,createWS,updateWS!,startingPointWS,optimise,residual,Fgradient,hessianApprox,points,descent,levMar,optim,lbfgsb,wise,readFilter,readWS

const R = Float64;
const T = Complex{R};

mutable struct WeightSearch
    intv::Vector{R}
    intv_free
    intv_free_count
    wght::Vector{R}
    wght_free
    wght_free_count
    intv_modifier
end

function createWS( intv,
                   intv_free::BitArray{1},
                   wght,
                   wght_free::BitArray{1},
                   intv_modifier=x->x )
    # Check vector sizes

    if !(length(intv) == length(intv_free) && length(wght) == length(wght_free) && length(intv) == length(wght))
        error("Invalid vector input.\nIt must hold that length(intv) == length(intv_free) && length(wght) == length(wght_free) && length(intv) == length(wght).")
    end

    return SLiSeFilters.WeightSearch(
        vec(intv),
        intv_free,
        length(intv[intv_free]),
        vec(wght),
        wght_free,
        length(wght[wght_free]),
        intv_modifier
    )
end

function startingPointWS( ws::SLiSeFilters.WeightSearch )
    return vcat( ws.intv[ws.intv_free],ws.wght[ws.wght_free] )
end

function updateWS!( ws::SLiSeFilters.WeightSearch, x::Vector )
    ws.intv[ws.intv_free] = x[1:ws.intv_free_count]
    ws.wght[ws.wght_free] = x[ws.intv_free_count+1:end]
    ws.intv_modifier(ws.intv)
    return SLiSeFilters.createConfig( ws.intv, ws.wght );
end

mutable struct Config
    Weights::Array{R,2}    # contains an interval a,b and the corresponding weight w in each column [a;b;w]
    Weights_h::Array{R,2}  # contains the sam as above e for the step-function h --> only in the interval [-1,1]
end

# Some predefined weights
function createConfig( name::String = "gamma" )
    # Gamma-SLISE as per paper
    if (name == "gamma")
        intgl = [ 0.95,  1.05, 1.4, 10];
        wght =  [1, .01, 10, 20];

        # Zeta-SLISE
    elseif (name == "zeta")
        intgl = [   0.95,    0.995,    1.005,     1.05,     1.1,   1.2,    5];
        wght  = [1,       4,        0.5,     4,        1.7,      1,   0.35];

        # Box-SLISE as per paper
    elseif (name == "box")
        intgl = [   0.95,    0.995,    1.005,     1.05,     1.1,   1.3, 1.8,   3];
        wght  = [      1,        4,        2,         4,    0.6,     1, 0.3, 0.1];
    end

    return createConfig(intgl,wght)
end

function createConfig( intgl::Array{R,1}, wght::Array{R,1} )
    intgl_h = intgl;
    wght_h = wght;

    intgl = vcat(-reverse(intgl), intgl);
    wght = vcat(reverse(wght), wght[2:end]);
    intervals = [intgl[1:end-1] intgl[2:end] wght]';

    intgl_h = vcat(intgl_h[ intgl_h .< 1.0], 1);
    wght_h = wght_h[ 1:length( intgl_h ) ];
    intgl_h = vcat(-reverse(intgl_h), intgl_h);
    wght_h = vcat(reverse(wght_h), wght_h[2:end]);
    intervals_h = [intgl_h[1:end-1] intgl_h[2:end] wght_h]';

    return Config(
        intervals,
        intervals_h
    );

end

function updateConfig!( c::SLiSeFilters.Config, intgl::Array{R,1}, wght::Array{R,1} )
    intgl_h = intgl;
    wght_h = wght;

    intgl = vcat(-reverse(intgl), intgl);
    wght = vcat(reverse(wght), wght[2:end]);
    c.Weights = [intgl[1:end-1] intgl[2:end] wght]';

    intgl_h = vcat(intgl_h[ intgl_h .< 1.0], 1);
    wght_h = wght_h[ 1:length( intgl_h ) ];
    intgl_h = vcat(-reverse(intgl_h), intgl_h);
    wght_h = vcat(reverse(wght_h), wght_h[2:end]);
    c.Weights_h = [intgl_h[1:end-1] intgl_h[2:end] wght_h]';

end

function optimise( methodName::String, z, alpha, c, lb = nothing )
    # Cast string to function
    f = getfield(SLiSeFilters, Symbol(methodName))

    if (lb == nothing)
        return f( z,alpha,c )
    else
        try
            return f( z,alpha,c,lb )
        catch e
            println(e)
            return lbfgsb( z,alpha,c,lb )
        end
    end
end

include("./definiteIntegrals.jl")
include("./gradMatrices.jl")
include("./residual.jl")
include("./points.jl")

include("./descent.jl")
include("./levMar.jl")
include("./lbfgsb.jl")
include("./optim.jl")

include("./wcr.jl")
include("./wise.jl")
include("./io.jl")

end # module
