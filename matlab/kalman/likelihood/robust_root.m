function RootV = robust_root( V )

    [ U, D ] = schur( 0.5 * ( V + V.' ) );
    d = diag( D );
    assert( max( max( abs( D - diag( d ) ) ) ) < eps );
    assert( isreal( U ) );
    
    Select = d > 0;
    
    RootV = U( :, Select ) * diag( sqrt( d( Select ) ) );

end
