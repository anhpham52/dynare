function hessian_mat = outer_product_gradient( func, x, dataset, varargin )

    if ~isa( func, 'function_handle' )
        func = str2func( func );
    end
    
    data = dataset.data;
    
    prior_weights = sum( ~isnan( data ), 2 );
    prior_weights = prior_weights ./ sum( prior_weights );

    dscores = GetJacobian( @( x ) get_lik_components( x, prior_weights, func, dataset, varargin ), x, length( prior_weights ) );
    
    hessian_mat = dscores' * dscores;

    hessian_mat = hessian_mat(:)';

end

function lik_components = get_lik_components( x, prior_weights, func, dataset, varargin )

    [ ~, ~, ~, lik_components ] = func( x, dataset, varargin{:} );
    
    assert( ~isempty( lik_components ) );
    prior = lik_components( 1 );
    lik_components( 1 ) = [];
    assert( length( lik_components ) == length( prior_weights ) );
    lik_components = lik_components(:) + prior_weights .* prior;

end
