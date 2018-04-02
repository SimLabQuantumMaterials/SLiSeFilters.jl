module SLiSeFilters
export residual,createConfig,Fgradient,descent,levMar,optimise,points

const R = Float64;
const T = Complex{R};

mutable struct Config
    # contains an interval a,b and the corresponding weight w in each column [a;b;w]
    Weights::Array{R,2}
    Weights_h::Array{R,2}
  # contains the same as above for the step-function h --> only in the interval [-1,1]
end

function createConfig(name::String="gamma")
    # Gamma-SLISE as per paper
    intgl = [ 0.95,  1.05, 1.4, 10];
    wght =  [1, .01, 10, 20];

    # OLD Gamma-SLISE
    if (name == "gamma_old")
        intgl = [ 0.95,  1.05, 1.4, 10];
        wght =  [50, 5.0, 500, 5000];

        # Zeta-SLISE
    elseif (name == "zeta")
        intgl = [   0.95,    0.995,    1.005,     1.05,     1.1,   1.2,    5];
        wght  = [1,       4,        0.5,     4,        1.7,      1,   0.35];

        # Box-SLISE as per paper
    elseif (name == "box")
        intgl = [   0.95,    0.995,    1.005,     1.05,     1.1,   1.3, 1.8,   3];
        wght  = [      1,        4,        2,         4,    0.6,     1, 0.3, 0.1];

        # E-Gamma-Slise
    elseif (name == "egamma")
        intgl = [0.96,  0.96^-1, 1.4, 10];
        wght  = [0.7, .00092, 887, 20];
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

include("./definiteIntegrals.jl")
include("./gradMatrices.jl")
include("./residual.jl")
include("./points.jl")
include("./descent.jl")
include("./levMar.jl")

end # module
