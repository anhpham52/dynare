function [ Residual, Jacobian ] = DynamicResidualWithGrowth( xRelLagSteadyNew, xRelLagSteady, xRelLagSteadyNewIndices, fname, dr, k2, iyr0, exo_simul, params, it_ )

    xRelLagSteady( xRelLagSteadyNewIndices ) = xRelLagSteadyNew;
    
    order_var = dr.order_var;

    xRelCurrentSteady              = zeros( size( xRelLagSteady ) );
    xRelCurrentSteady( order_var ) = dr.non_bgp_drift + dr.ghx * xRelLagSteady( order_var( k2 ) );
    xRelFutureSteady               = zeros( size( xRelLagSteady ) );
    xRelFutureSteady( order_var )  = dr.non_bgp_drift + dr.ghx * xRelCurrentSteady( order_var( k2 ) );
    
    z = [ xRelLagSteady + dr.ys; xRelCurrentSteady + dr.ys; xRelFutureSteady + dr.ys ];
    
    if nargout > 1
        
        [ Residual, jacobia_ ] = feval( [ fname '_dynamic' ], z( iyr0 ), exo_simul, params, dr.ys, it_ );

        d_xRelLagSteady_d_xRelLagSteadyNew     = zeros( length( xRelLagSteady ), length( xRelLagSteadyNew ) );
        d_xRelCurrentSteady_d_xRelLagSteadyNew = d_xRelLagSteady_d_xRelLagSteadyNew;
        d_xRelFutureSteady_d_xRelLagSteadyNew  = d_xRelLagSteady_d_xRelLagSteadyNew;

        d_xRelLagSteady_d_xRelLagSteadyNew( xRelLagSteadyNewIndices, : ) = eye( length( xRelLagSteadyNew ) );
        d_xRelCurrentSteady_d_xRelLagSteadyNew( order_var, : )        = dr.ghx * d_xRelLagSteady_d_xRelLagSteadyNew( order_var( k2 ), : );
        d_xRelFutureSteady_d_xRelLagSteadyNew( order_var, : )         = dr.ghx * d_xRelCurrentSteady_d_xRelLagSteadyNew( order_var( k2 ), : );

        Jacobian = zeros( length( Residual ), length( z ) );
        Jacobian( :, iyr0 ) = jacobia_( :, 1 : length( iyr0 ) );
        Jacobian = Jacobian * [ d_xRelLagSteady_d_xRelLagSteadyNew; d_xRelCurrentSteady_d_xRelLagSteadyNew; d_xRelFutureSteady_d_xRelLagSteadyNew ];
    
    else
        
        Residual = feval( [ fname '_dynamic' ], z( iyr0 ), exo_simul, params, dr.ys, it_ );

    end

end
