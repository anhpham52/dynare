function Jacobian = GetJacobian( f, x, nf, h )
    nx = length( x );
    Jacobian = NaN( nf, nx );
    sreps = sqrt( eps );
    if nargin < 4
        h = eps^(1/3);
    end
    h = ones( nx, 1 ) .* h;
    parfor i = 1 : nx
        WarningState = warning( 'off', 'all' );
        xi = x( i );
        hi = max( [ h( i ), sreps, abs( sreps * xi ) ] );
        while true
            BreakFlag = false;
            if hi < eps
                hi = eps;
                BreakFlag = true;
            end
            if hi < eps( xi )
                xi = eps( xi );
                BreakFlag = true;
            end
            try
                Jacobian( :, i ) = ( f( SetElement( x, i, xi + hi ) ) - f( SetElement( x, i, xi - hi ) ) ) / ( 2 * hi ); %#ok<PFBNS>
            catch Error
                DisplayError( Error );
            end
            if all( isfinite( Jacobian( :, i ) ) ) || BreakFlag
                break
            else
                hi = 0.5 * hi;
            end
        end
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
