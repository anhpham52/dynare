function [ info, M, options, oo ] = non_central_approximation( M, options, oo )
    options.original_order = options.order;
    [ deflect, M, oo ] = compute_deflected_linear_approximation( M, options, oo, options.non_central_approximation );
    options.order = 1;
    oo.dr.ys = deflect.y;
    oo.steady_state = deflect.y;
    oo.dr.ghx = deflect.y_x;
    oo.dr.ghu = deflect.y_u;
    
    info = 0;
    if any( isnan( deflect.y ) )
        info = 12345;
    end
    if any( any( isnan( deflect.y_x ) ) )
        info = 12345;
    end
    if any( any( isnan( deflect.y_u ) ) )
        info = 12345;
    end
end