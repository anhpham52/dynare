function out = logdet( in )
    % det( in ) = prod( eig( in ) )
    % log( det( in ) ) = sum( log( eig( in ) ) )
    out = sum( log( max( realmin, eig( in ) ) ) );
end
