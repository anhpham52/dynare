function h = FindStepSize( f, x, h )
    nx = length( x );
    sreps = sqrt( eps );
    if nargin < 3
        h = eps^(1/3);
    end
    h = ones( nx, 1 ) .* h;
    fx = f( x );
    parfor i = 1 : nx
        WarningState = warning( 'off', 'all' );
        xi = x( i );
        hi = max( [ h( i ), sreps, abs( sreps * xi ), eps( abs( xi ) ) ] );
        while true
            if ( f( SetElement( x, i, xi + hi ) ) > fx ) && ( f( SetElement( x, i, xi - hi ) ) > fx ) %#ok<PFBNS>
                break
            end
            hi = 2 * hi;
        end
        h( i ) = hi;
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
