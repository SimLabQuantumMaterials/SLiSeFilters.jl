import Optim

# Transformation from C to R
split_(x)  = x[1:Int(length(x)/2)],x[Int(length(x)/2) + 1:end]

"BFGS using Hager-Zhang"
function optim( z, alpha, c, lb=nothing; linesearch = Optim.LineSearches.MoreThuente() )
    function grad!(stor,x)
        stor[:] = vcat(Fgradient(split_(x)..., c)...)
    end

    result = Optim.optimize( x->residual(split_(x)...,c),
                             grad!,
                             vcat(z, alpha),
                             Optim.BFGS(linesearch = linesearch), Optim.Options(
                                iterations = 14000,
                                allow_f_increases = true))

    return split_(Optim.minimizer(result))
end
