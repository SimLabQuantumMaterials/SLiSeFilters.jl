"Evaluates a rational filter, exploiting numerical optimizations provided by JULIA."
function phi(x::Float64, z, alpha)
    y = zero(Complex{Float64});
    @inbounds @fastmath for i = eachindex(z) # @fastmath can be dangerous but no impact has been noticeable in various tests
        y = y + alpha[i] / (x[1] - z[i]) + conj(alpha[i]) / (x[1] - conj(z[i])) - alpha[i] / (x[1] + z[i]) - conj(alpha[i]) / (x[1] + conj(z[i]));
    end
    return real(y)
end

# For compatibility with JULIA versions => 0.7
function Linspace(left, right, points)
    return left:(right-left)/(points-1):right
end


"""
    wcr( z, a, G; interval_fevals=500000, upper=10., middle=2.)

Obtain the WCR of a rational filter.

Samples thousands points, exploiting numerical optimizations provided by JULIA
in order to obtain the WCR.

In terms of both computational time and, most importantly, validity, the employed
approach is the only feasible choice, as opposed to sampling the vector-based phi()
function or using the eigenvalues of the companion matrix of the nominator of the
rational filter or of its derivative, whilst both Float64 and BigFloat were tested
but could not overcome the occuring large condition numbers. This has not produced
any large derivations from the expected WCR on a set of 60 'difficult' rational
filters. More such filters are to be obtained throghout the subsequent tests in
filter computation.

Underlying assumptions: 1. Rational filters are continuous, thus limited in terms
of derivation. 2. We study at worst G=0.99998, thus we want to be able to consider
lenghts less than 0.99998^-1-1, which requires at least (0.99998^-1-1)^-1 = 49999
evaluations in (1,2)."""
function wcr(z,a,G; interval_fevals=500000, upper=10., middle=2., lowerGap = G, mintol = Inf)
    # Compute maximizer outside (-1,1) with respect to the gap
    # Proceed in three steps, contributing different precision to each interval
    maximizer = 0
    for x in Linspace(G^-1,1.1+(G^-1-1),interval_fevals)
        w = abs(phi(x,z,a))
        if (w > maximizer)
            maximizer = w
        end
    end
    for x in Linspace(1.1+(G^-1-1),middle+(G^-1-1),interval_fevals)
        w = abs(phi(x,z,a))
        if (w > maximizer)
            maximizer = w
        end
    end
    for x in Linspace(middle+(G^-1-1),upper+(G^-1-1),interval_fevals)
        w = abs(phi(x,z,a))
        if (w > maximizer)
            maximizer = w
        end
    end

    # Compute the minimizer inside (-1,1)
    minimizer = mintol
    for x in Linspace(0,lowerGap,interval_fevals)
        w = abs(phi(x,z,a))
        if (w < minimizer)
            minimizer = w
        end
    end

    return maximizer/minimizer
end
