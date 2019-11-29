import Lbfgsb

# 3<=m<=20 recommended
# factr =1.d+12 for low accuracy; factr =1.d+7 for moderate accuracy; factr =1.d+1 for extremely high accuracy.
# If factr =0, the test will stop the algorithm only if the objective function remains unchanged after one iteration.
function lbfgsb( z, alpha, c, lb = nothing; factr = 1e1, m=30, maxiter = 5000, iprint=-1, pgtol=1e-7)
  x       = vcat(z, alpha)
  x_      = [real(x);imag(x)]
  n       = length(x_)
  
  function objFunc(x)
    z, alpha = unpackComplex(x)
    z = residual(z, alpha, c)
    return z
  end
  
  function gradFunc!(x, g::Array)
    z, alpha = unpackComplex(x)
    y = vcat(Fgradient(z,alpha, c)...)
    g[:] = 2*[real(y);imag(y)]
    return
  end
  
  if (lb != nothing)
    lower = [ones(div(n,2))*-Inf;ones(div(n,4))*lb;ones(div(n,4))*-Inf] # Stores the lower bounds for every entry of x
    btype = [ones(div(n,2))*0;ones(div(n,4))*1;ones(div(n,4))*0] # Set bounds
    x_[x_ .< lower] = lower[x_ .< lower] # project into bounds
    if (!isempty(x_[x_ .< lower]))
      error("Could not project into boundaries.")
    end
    f, x, numCall, numIter, status = Lbfgsb.lbfgsb( objFunc, gradFunc!, x_, btype=btype, lb=lower, factr=factr, m=m, maxiter=maxiter, iprint=iprint, pgtol=pgtol )
  else
    f, x, numCall, numIter, status = Lbfgsb.lbfgsb( objFunc, gradFunc!, x_, m=m, factr=factr, maxiter=maxiter, iprint=iprint, pgtol=pgtol )
  end
  
  return unpackComplex(x)
end

function splitVec( x )
  return x[1:Int(length(x)/2)],x[Int(length(x)/2) + 1:end]
end

function unpackComplex(x)
  Re,Im = splitVec(x)
  z_re,a_re = splitVec(Re)
  z_im,a_im = splitVec(Im)
  
  return z_re + im .* z_im, a_re + im .* a_im
end
