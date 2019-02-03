function HessianDiagonal = GetHessianDiagonal( f, x, nf, EnsureSign, h )
    nx = length( x );
    HessianDiagonal = NaN( nf, nx );
    sreps = sqrt( eps );
    if nargin < 5
        h = eps^(1/6);
    end
    h = ones( nx, 1 ) .* h;
    fx = f( x );
    assert( all( isfinite( fx ) ) );
    parfor i = 1 : nx
        fprintf( '\nStarting obtaining second derivative of coordinate %d.\n', i );
        WarningState = warning( 'off', 'all' );
        xi = x( i );
        hi = max( [ h( i ), sreps, abs( sreps * xi ) ] );
        Expanded = false;
        Shrank = false;
        while true
            if hi < sreps
                hi = sreps;
                Expanded = true;
            end
            if hi < sqrt( eps( xi ) )
                hi = sqrt( eps( xi ) );
                Expanded = true;
            end
            if hi >= max( abs( xi ), 1 )
                hi = max( abs( xi ), 1 );
                Shrank = true;
            end
            fp2 = NaN( size( fx ) );
            fp  = NaN( size( fx ) );
            fn  = NaN( size( fx ) );
            fn2 = NaN( size( fx ) );
            try
                fp = f( SetElement( x, i, xi + hi ) ); %#ok<PFBNS>
                if all( isfinite( fp ) )
                    try
                        fp2 = f( SetElement( x, i, xi + 2 * hi ) );
                    catch Error
                        DisplayError( Error );
                    end
                end
            catch Error
                DisplayError( Error );
            end
            try
                fn = f( SetElement( x, i, xi - hi ) );
                if all( isfinite( fn ) )
                    try
                        fn2 = f( SetElement( x, i, xi - 2 * hi ) );
                    catch Error
                        DisplayError( Error );
                    end
                end
            catch Error
                DisplayError( Error );
            end
            d = zeros( nf, 5 );
            d( :, 1 ) = ( -fn2 / 12 + 4 / 3 * fn - 5 / 2 * fx + 4 / 3 * fp - fp2 / 12 ) / ( hi * hi );
            d( :, 2 ) = ( fn  - 2 * fx + fp  ) / ( hi * hi );
            d( :, 3 ) = ( fn2 - 2 * fx + fp2 ) / ( 4 * hi * hi );
            d( :, 4 ) = ( fx  - 2 * fp + fp2 ) / ( hi * hi );
            d( :, 5 ) = ( fn2 - 2 * fn + fx  ) / ( hi * hi );
            Possible = all( isfinite( d ), 1 );
            if EnsureSign == 0
                Preferable = true( 1, 5 );
            else
                Preferable = all( sign( d ) == EnsureSign, 1 );
            end
            Good = find( Possible & Preferable, 1 );
            OK   = find( Possible, 1 );
            if ~isempty( Good )
                HessianDiagonal( :, i ) = d( :, Good );
                fprintf( '\nCompleted obtaining second derivative of coordinate %d.\n', i );
                break
            elseif ~isempty( OK )
                HessianDiagonal( :, i ) = d( :, OK );
                if Shrank
                    break
                end
                hi = 4 * hi;
                Expanded = true;
                fprintf( '\nExpanding coordinate %d.\n', i );
            else
                if Expanded
                    break
                end
                hi = 0.25 * hi;
                Shrank = true;
                fprintf( '\nShrinking coordinate %d.\n', i );
            end
        end
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
