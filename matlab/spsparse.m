function x = spsparse( x )

    x = squeeze( x );
    
    sx = size( x );
    
    if length( sx ) > 2
        
        dimDist = arrayfun( @( d ) ones( d, 1 ), sx( 3 : end ), 'UniformOutput', false );
        
        x = squeeze( mat2cell( x, sx( 1 ), sx( 2 ), dimDist{:} ) );
        
        x = cellfun( @( x_ ) spsparse( x_ ), x, 'UniformOutput', false );
        
    else

        x = sparse( x );
        x( abs( x ) < 1.81898940354586e-12 ) = 0;
    
    end

end
