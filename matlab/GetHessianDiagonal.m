function HessianDiagonal = GetHessianDiagonal( f, x, nf, EnsureSign, h )
    nx = length( x );
    HessianDiagonal = NaN( nf, nx );
    sreps = sqrt( eps );
    if nargin < 5
        h = eps^(1/3);
    end
    h = ones( nx, 1 ) .* h;
    fx = f( x );
    assert( all( isfinite( fx ) ) );
    parfor i = 1 : nx
        WarningState = warning( 'off', 'all' );
        xi = x( i );
        hi = max( [ h( i ), sreps, abs( sreps * xi ) ] );
        Expanded = false;
        Shrank = false;
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
                if all( isfinite( fp ) )
                    try
                        fp2 = f( SetElement( x, i, xi + 2 * hi ) );
                    catch Error
                        DisplayError( Error );
                        fp2 = NaN( size( fx ) );
                    end
                else
                    fp2 = NaN( size( fx ) );
                end
            catch Error
                DisplayError( Error );
                fp = NaN( size( fx ) );
                fp2 = NaN( size( fx ) );
            end
            try
                fn = f( SetElement( x, i, xi - hi ) );
                if all( isfinite( fn ) )
                    try
                        fn2 = f( SetElement( x, i, xi - 2 * hi ) );
                    catch Error
                        DisplayError( Error );
                        fn2 = NaN( size( fx ) );
                    end
                else
                    fn2 = NaN( size( fx ) );
                end
            catch Error
                DisplayError( Error );
                fn = NaN( size( fx ) );
                fn2 = NaN( size( fx ) );
            end
            GoodFlag = true;
            if all( isfinite( fp ) )
                if all( isfinite( fp2 ) )
                    if all( isfinite( fn ) )
                        if all( isfinite( fn2 ) )
                            HessianDiagonal( :, i ) = ( -fn2 / 12 + 4 / 3 * fn - 5 / 2 * fx + 4 / 3 * fp - fp2 / 12 ) / ( hi * hi );
                        else
                            HessianDiagonal( :, i ) = ( fn - 2 * fx + fp ) / ( hi * hi );
                        end
                    else
                        HessianDiagonal( :, i ) = ( fx - 2 * fp + fp2 ) / ( hi * hi );
                    end
                else
                    if all( isfinite( fn ) )
                        HessianDiagonal( :, i ) = ( fn - 2 * fx + fp ) / ( hi * hi );
                    else
                        GoodFlag = false;
                    end
                end
            else
                if all( isfinite( fn ) )
                    if all( isfinite( fn2 ) )
                        HessianDiagonal( :, i ) = ( fn2 - 2 * fn + fx ) / ( hi * hi );
                    else
                        GoodFlag = false;
                    end
                else
                    GoodFlag = false;
                end
            end
            if GoodFlag 
                if ( EnsureSign == 0 ) || all( sign( HessianDiagonal( :, i ) ) == EnsureSign )
                    break
                else
                    if Shrank
                        break
                    end
                    hi = 2 * hi;
                    Expanded = true;
                    if hi > max( abs( xi ), 1 )
                        break
                    end
                end
            elseif BreakFlag
                break
            else
                if Expanded
                    break
                end
                hi = 0.5 * hi;
                Shrank = true;
            end
        end
        warning( WarningState );
    end
end

function x = SetElement( x, i, xi )
    x( i ) = xi;
end
