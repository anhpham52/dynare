function Jacobian = GetJacobian( f, x, nf, h )
    nx = length( x );
    Jacobian = NaN( nf, nx );
    sreps = sqrt( eps );
    if nargin < 4
        h = eps^(1/3);
    end
    h = ones( nx, 1 ) .* h;
    fx = f( x );
    assert( all( isfinite( fx ) ) );
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
                fp = f( SetElement( x, i, xi + hi ) ); %#ok<PFBNS>
            catch Error
                DisplayError( Error );
                fp = NaN( size( fx ) );
            end
            try
                fn = f( SetElement( x, i, xi - hi ) );
            catch Error
                DisplayError( Error );
                fn = NaN( size( fx ) );
            end
            if all( isfinite( fp ) )
                if all( isfinite( fn ) )
                    Jacobian( :, i ) = ( fp - fn ) / ( 2 * hi );
                    break
                else
                    Jacobian( :, i ) = ( fp - fx ) / hi;
                    break
                end
            else
                if all( isfinite( fn ) )
                    Jacobian( :, i ) = ( fx - fn ) / hi;
                    break
                else
                    hi = 0.5 * hi;
                end
            end
            if BreakFlag
                break
            end
        end
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
