import BlackBoxOptim
import Optim
import NLopt

"Perform WCR minimization, through the SLiSe framework"
function wise( ws, z1, a1, G_user;
            G_in = 1,
            G_out = 2,
            outpath = nothing,
            lb = nothing,
            refine = true,

            # To disable parts of oriented restart
            gappingEnabled = true,
            weightsEnabled = true,
            lastIntEnabled = true
        )
    p      = length( z1 )
    x0     = SLiSeFilters.startingPointWS(ws) # the starting weight function
    G      = sqrt( G_user )                   # shifting for better filters
    mintol = Inf

    # Use bounds within minimization
    minalgo = "optim"
    if lb != nothing
        minalgo = "lbfgsb"
    end

    # Maximum iterations of minimizations
    maxiter          = 200
    every            = 10
    runs             = 20
    inner_max_eval   = 300
    inner_max_eval_2 = 1000

    # Strech variables for parameters in blackbox optimization
    weight_multiplier    = 10.
    interval_multplier   = 3.

    # Remeber set-up, to be able to change and reset it easily
    intv_free_prev       = copy(ws.intv_free)
    intv_free_count_prev = copy(ws.intv_free_count)
    wght_free_prev       = copy(ws.wght_free)
    wght_free_count_prev = copy(ws.wght_free_count)

    # Monitor progress
    prevWCR    = Inf # monitor improvement step by step
    prevItWcr  = Inf # monitor improvement iteration by iteration
    tol        = 1e-9
    tolReached = false

    # Objective function
    function f(value::Vector)
        c      = SLiSeFilters.updateWS!(ws, abs.(value))
        z3, a3 = SLiSeFilters.optimise( minalgo, z1, a1, c, lb )

        # flush output to console
        if VERSION >= v"0.7.0"
            flush(stdout)
        else
            flush(STDOUT)
        end

        return wcr( z3, a3, G, mintol = mintol, lowerGap = G );
    end

    # For oriented restart
    function callBlackBox(ws, z1, a1, prevWCR, searchRange, i; traceMode=:compact, evals=inner_max_eval, method=:adaptive_de_rand_1_bin_radiuslimited) # method is default method in BlackBoxOptim
        x0  = SLiSeFilters.startingPointWS(ws)
        res = BlackBoxOptim.bboptimize(f; SearchRange = searchRange, NumDimensions = length(searchRange), MaxFuncEvals = evals, Method=method, TraceMode=traceMode, TraceInterval=20.0);
        bf  = BlackBoxOptim.best_fitness(res)
        bs  = BlackBoxOptim.best_candidate(res)

        # Paramter Gap: Evaluation and reset
        if (bf < prevWCR)
            @printf "[%d] Setting to new min: %.16e\n" i bf
            prevWCR = bf

            c       = SLiSeFilters.updateWS!( ws, abs.(bs) )
            z1,a1   = SLiSeFilters.optimise( minalgo, z1, a1, c, lb )
        else
            SLiSeFilters.updateWS!(ws, x0) # to call modifier before resetting it
        end

        return z1,a1,prevWCR
    end

    for i = 1:runs
        println("Starting.")

        ws.wght_free       = falses(length(ws.wght_free))
        ws.intv_free       = falses(length(ws.intv_free))
        ws.intv_free_count = 0
        ws.wght_free_count = 0

        # Oriented restart/Minimize gap position
        if (gappingEnabled)
            ws.intv_free_count = 1

            @printf "Performing gap minimization [1/2]..\n"
            ws.intv_free[G_in] = true
            z1,a1,prevWCR = callBlackBox(ws, z1, a1, prevWCR, [(G_user, 1.)], i)
            ws.intv_free[G_in]  = false

            @printf "Performing gap minimization [2/2]..\n"
            ws.intv_free[G_out] = true
            z1,a1,prevWCR = callBlackBox(ws, z1, a1, prevWCR, [(1., G_user^-1)], i)
            ws.intv_free[G_out] = false

            ws.intv_free_count = 0
        end

        # Oriented restart/Loop through weights (except gaps)
        if (weightsEnabled)
            ws.wght_free_count  = 1

            for j = 1:length(ws.wght_free)
                @printf "Performing weight (%d) minimization..\n" j
                ws.wght_free[j] = true
                x0  = SLiSeFilters.startingPointWS(ws)
                z1,a1,prevWCR = callBlackBox(ws, z1, a1, prevWCR, [(x0[1]/weight_multiplier, weight_multiplier*x0[1])], i)
                ws.wght_free[j] = false
            end

            ws.wght_free[end]  = false
            ws.wght_free_count = 0
        end

        # Oriented restart/Special, dependent minimization for last interval
        if (lastIntEnabled)
            ws.wght_free_count = 1
            ws.intv_free_count = 1
            ws.wght_free[end]  = true
            ws.intv_free[end]  = true

            @printf "Performing last interval minimization..\n"
            x0  = SLiSeFilters.startingPointWS(ws)
            z1,a1,prevWCR = callBlackBox(ws, z1, a1, prevWCR, [(ws.intv[end-1], interval_multplier*ws.intv[end]), (x0[2]/weight_multiplier, weight_multiplier*x0[2])], i, evals=inner_max_eval_2, method=:resampling_memetic_search)
        end

        # Reset free variables
        ws.intv_free       = intv_free_prev
        ws.wght_free       = wght_free_prev
        ws.intv_free_count = intv_free_count_prev
        ws.wght_free_count = wght_free_count_prev

        # Do Nelder-Mead
        x0  = SLiSeFilters.startingPointWS(ws)
        res = Optim.optimize(f,
	                         x0,
	                         Optim.NelderMead( ),
	                         Optim.Options(iterations        = maxiter,
	                                       g_tol             = 1e-10, # =√(Σ(yᵢ-ȳ)²)/n, the conv crit for NelderMead
	                                       show_every        = every,
	                                       allow_f_increases = true, # in the end, minimizer is returned
	                                       show_trace        = true
	                                      )
	                        )

        # Evaluation and reset
        if (res.minimum < prevWCR)
            @printf "[%d] Setting to new min: %.16e\n" i res.minimum
            prevWCR = res.minimum

            c       = SLiSeFilters.updateWS!( ws, abs.(res.minimizer) )
            z1,a1   = SLiSeFilters.optimise( minalgo, z1, a1, c, lb )
        else
            SLiSeFilters.updateWS!(ws, x0) # reset search space
        end

        # Cache the results
        if (outpath != nothing)
            writedlm( string(outpath,"/wcr_",p,"_",G_user,"_lb",lb,"_filter_it",i,"_wise.dat"),([z1 a1]/G)',',') # comma-separated, for JULIA v0.7
            open( string(outpath,"/wcr_",p,"_",G_user,"_lb",lb,"_filter_it",i,"_wise.dat"), "a") do f # append
                writecsv(f, ([real(z1) imag(z1) real(a1) imag(a1)]/G)') # comma-separated, for JULIA v0.6
            end
            writedlm(string(outpath,"/wcr_",p,"_",G_user,"_lb",lb,"_weight_it",i,"_wise.dat"),[ws.intv ws.wght]',',')
        end

        # Validate change in WCR
        resd = abs( prevItWcr - prevWCR )/prevWCR
        @printf "[%d] CurrentMin: %.16e (internal)\t%.16e (real)" i prevWCR wcr( z1/G, a1/G, G_user ) # It holds that wcr( z1, a1, G, lowerGap = G ) = wcr( z1/G, a1/G, G_user, lowerGap = 1. )
        if (resd < tol)
            # If no progress in two subsequent runs, finish
            if (tolReached)
                @printf "\tEND. TOL reached.\n"
                break
            else
                @printf "\tTOL reached. Continuing..\n"
                tolReached = true
            end
        else
            tolReached = false
            if (i == runs)
                @printf "\tEND. MAXITER reached."
            end
            @printf "\n"
        end
        prevItWcr = prevWCR
    end

    if (!refine)
        return z1, a1
    end

    @printf "Performing refinement step..\n"
    return wiseRefine( z1, a1, G_user, outpath = outpath )
end

"Perform WCR minimization, without SLiSe"
function wiseRefine( z1, a1, G_user; lb = nothing, outpath = nothing )
    p       = length( z1 )
    maxiter = 100000
    every   = 100
    mintol  = Inf

    # Transformation from C to R
    toR(x)     = [real(x); imag(x)]
    fromR(x)   = x[1:Int(length(x)/2)] + im*x[Int(length(x)/2) + 1:end]
    x0         = toR([z1; a1])

    function g(x)
        # flush output to console
        if VERSION >= v"0.7.0"
            flush(stdout)
        else
            flush(STDOUT)
        end
        return wcr(split_(fromR(x))...,G_user, mintol = mintol )
    end

    count = 0
    function f2(value::Vector, grad::Vector)
        count::Int += 1

        y = wcr(split_(fromR(value))...,G_user, mintol = mintol )

        if (count % every == 0)
         @printf "%6d   %14e\n" count y
        end

        return y
    end

    # Do the optimization
    z1, a1 = nothing, nothing

    if (lb != nothing)
        n     = length(x0)
        lower = [ones(div(n,2))*-Inf;ones(div(n,4))*lb;ones(div(n,4))*-Inf]
        upper = [ones(n) * Inf]

        opt = NLopt.Opt(:LN_NELDERMEAD, length(x0))
        NLopt.lower_bounds!( opt, lower )
        NLopt.maxeval!(opt, maxiter)
        NLopt.ftol_rel!(opt, 1e-6)
        NLopt.min_objective!(opt, f2)

        # Perform NLOpt Minimization
        @printf "Evals    Function value    \n"
        @printf "------   --------------    \n"
        (minf,minx,reason) = NLopt.optimize(opt, x0)

        z1,a1 = split_(fromR(minx))
    else
        res = Optim.optimize(g,
                   x0,
                   Optim.NelderMead(),
                   Optim.Options(iterations        = maxiter,
                                 g_tol             = 1e-13, # =√(Σ(yᵢ-ȳ)²)/n, the conv crit for NelderMead
                                 show_every        = every,
                                 allow_f_increases = true, # in the end, minimizer is returned
                                 show_trace        = true
                                 )
                   )

        # Check for convergence
        if !Optim.converged(res)
            @printf "Attention: No convergence reached after %d steps!\n" maxiter
        end

        # Obtain the results
        z1,a1 = split_(fromR(res.minimizer))
    end

    # Output the results
    @printf "[Refined] CurrentMin: %.16e\n" wcr( z1, a1, G_user )

    # Save the results to disk
    if (outpath != nothing)
        writedlm( string(outpath,"/wcr_",p,"_",G_user,"_filter_wise_refined.dat"),([z1 a1])',',') # comma-separated, for JULIA v0.7
        open( string(outpath,"/wcr_",p,"_",G_user,"_filter_wise_refined.dat"), "a") do f # append
            writecsv(f, ([real(z1) imag(z1) real(a1) imag(a1)])') # comma-separated, for JULIA v0.6
        end
    end

    return z1,a1
end
