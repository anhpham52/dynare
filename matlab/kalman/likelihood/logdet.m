function out = logdet( in )
    % det( in ) = prod( eig( in ) )
    % log( det( in ) ) = sum( log( eig( in ) ) )
    % disp( sum( eig( in ) < realmin ) * log( realmin ) );
    out = sum( log( eig( in ) ) );
end
