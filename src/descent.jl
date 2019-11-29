if VERSION >= v"0.7.0"
    using Printf
    using LinearAlgebra
end

function descent( z::Array{T,1}, alpha::Array{T,1}, c::SLiSeFilters.Config, lb::R=R(0.0) )
    
    project = function( z::Array{T,1} )
        mask = abs.(imag(z)) .< lb
        z[mask] = real(z[mask]) + lb*im*sign.(imag(z[mask]));
        z
    end
    
    maxiter = 20000;
    i = 1;
    
    while i < maxiter
        dz,da = Fgradient( z, alpha, c );
        
        a = 0.5;
        b = 0.5;
        t::Float64 = 1*b;
        
        resd = residual(z, alpha, c);
        #emm = dot(-([da; dz]), [da; dz]);
        emm = -norm([da; dz])^2;
        
        while( residual( project(z-t .* dz), alpha - t.*da, c) > real(resd + a * t * emm) )
            t = t*b;
            if 1e-16 > t
                println("done...")
                i=maxiter;
                break;
            end
        end
        
        alpha -= t*da;
        z = project(z - t*dz);
        
        
        if mod(i,500) == 0
            @printf "%.16e\t%d\n" resd i
        end
        i+=1
    end
    
    return z,alpha
end
