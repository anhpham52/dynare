function A = qr0( A, ~ )

    if issparse( A )
        A = qr( A, 0 );
    else
        A = triu( qr( A, 0 ) );
        [ m, n ] = size( A );
        if m > n
            A = A( 1 : n, : );
        end
    end

end
