# Do the integration of the approximation of h
# --> applies a function to the intervals in c.Weights
function case( f::Function, c::SLiSeFilters.Config )
    R::T = 0.0;
    intervals = c.Weights;

    for i = 1:size(intervals,2)
        a = intervals[1,i];
        b = intervals[2,i];
        w = intervals[3,i];

        R = R + w * f(a,b)
    end
    R
end

# Do the integration of h
# --> applies a function to the intervals in c.Weights_h
function case_h( f::Function, c::SLiSeFilters.Config )
    R::T = 0;
    intervals = c.Weights_h;

    for i = 1:size(intervals,2)
        a = intervals[1,i];
        b = intervals[2,i];
        w = intervals[3,i];

        R = R + w * f(a,b)
    end
    R
end


# 1 / (t - z)
function case1( z1::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        log( (b-z1) ./ (a-z1) )
    end
    return case_h( f, c::SLiSeFilters.Config )
end


# 1 / (t - z)^2
function case2( z1::T, c::SLiSeFilters.Config )
    f = function ( a::R, b::R )
        return (b-a) ./ ( (b-z1) .* (a-z1) )
    end
    return case( f, c::SLiSeFilters.Config )
end


# h(t) / (t - z)^2
function case2_h( z1::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        (b-a) ./ ( (b-z1) .* (a-z1) )
    end
    return case_h( f, c::SLiSeFilters.Config )
end


function case3( z1::T, z2::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        1 ./ ( z1-z2 ) .*
        (
         log( (b - z1 ) ./ (a - z1) ) +
         log( (a - z2 ) ./ (b - z2) )
         )
    end
    return case( f, c::SLiSeFilters.Config )
end


# 1 / (t-z1)^2 / (t-z2)
function case4( z1::T, z2::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        l1 = (b-z1) / (a-z1);
        l2 = (b-z2) / (a-z2);

        (b-a) / ( (z1-z2)*(b-z1)*(a-z1) ) - 1 / (z1-z2)^2 * ( log( l1 ) ) + 1 / (z1-z2)^2 * ( log( l2 ) ) 
    end
    return case( f, c::SLiSeFilters.Config )
end


# 1 / (t-z1)^3
function case5( z1::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        0.5* (b-a) * (b+a-2*z1) / ((b-z1)^2 * (a-z1)^2)
    end
    return case( f, c::SLiSeFilters.Config )
end


function case6( z1::T, z2::T, c::SLiSeFilters.Config )

    f = function ( a, b )
        (a-b)/(z1-z2)^2 * (  1/(a-z1)/(-b+z1) + 1 / (a-z2)/(-b+z2) ) + 2/(z1-z2)^3   * (  log( (a-z1)/(b-z1) ) + log( (b-z2)/(a-z2)  )  )
    end
    return case( f, c::SLiSeFilters.Config )
end


# 1/(t-z1)^4
function case7( z1::T, c::SLiSeFilters.Config )
    f = function ( a, b )
        1/3 * ( 1/(a-z1)^3 + 1/(-b+z1)^3 )
    end

    return case( f, c::SLiSeFilters.Config )
end


#  1/(t-z1)^4
function caseh_h( c::SLiSeFilters.Config )
    f = function ( a, b )
        0.5*abs(b-a);
    end
    return case_h( f, c::SLiSeFilters.Config )
end
