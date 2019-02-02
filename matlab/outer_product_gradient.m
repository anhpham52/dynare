function hessian_mat = outer_product_gradient( func, x, dataset, varargin )

    if ~isa( func, 'function_handle' )
        func = str2func( func );
    end
    
    data = dataset.data;
    
    prior_weights = sum( ~isnan( data ), 2 );
    prior_weights = prior_weights ./ sum( prior_weights );
    
    seps = sqrt( eps );
    sseps = sqrt( seps );

    dscores = GetJacobian( @( x ) get_lik_components( x, prior_weights, func, dataset, varargin{:} ), x, length( prior_weights ), max( sseps, sseps * abs( x ) ) );
    
    bad_params = find( ~( all( isfinite( dscores ) ) ) );
    
    if ~isempty( bad_params )
        disp( 'Results may be unreliable due to non-finite Hessian elements with respect to the following parameters:' );
        disp( bad_params );
        dscores( :, bad_params ) = 0;
    end
    
    hessian_mat = dscores' * dscores;

    hessian_mat = 0.5 * ( hessian_mat + hessian_mat' );
    
    hessian_mat = hessian_mat + seps * eye( size( hessian_mat ) );

    if all( isfinite( hessian_mat(:) ) )
        eig_hessian_mat = eig( hessian_mat );
        disp( 'Negative elements of the eigenvalues of the Hessian:' );
        disp( eig_hessian_mat( eig_hessian_mat < 0 ) );
        hessian_mat = NearestSPD( hessian_mat );
        eig_hessian_mat = eig( hessian_mat );
        disp( 'Negative elements of the eigenvalues of the modified Hessian:' );
        disp( eig_hessian_mat( eig_hessian_mat < 0 ) );
    end

    hessian_mat = hessian_mat(:)';

end

function lik_components = get_lik_components( x, prior_weights, func, dataset, varargin )

    [ ~, ~, ~, lik_components ] = func( x, dataset, varargin{:} );
    
    assert( ~isempty( lik_components ) );
    lik_components = lik_components(:);
    prior = lik_components( 1 );
    lik_components( 1 ) = [];
    assert( length( lik_components ) == length( prior_weights ) );
    lik_components = lik_components + prior_weights .* prior;

end
