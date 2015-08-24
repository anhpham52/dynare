function RV = parallel_wrapper( objective_function, XV, varargin )
    n = size( XV, 2 );
    RV = zeros( 1, n );
    parfor i = 1 : n
        R = objective_function( XV( :, i ), varargin{:} ); %#ok<PFBNS>
        if isempty( R )
            R = 1e12;
        end
        RV( i ) = R;
    end
end
