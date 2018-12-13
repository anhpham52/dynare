function ParallelInitializationInternal( varargin )
    Seed = randi( ceil( intmax / 2 ) );
    if nargin == 0
        spmd
            SetGlobalWorkerNumber( labindex );
            rng( Seed + labindex );
        end
    else
        in_options_ = varargin{ 1 };
        spmd
            SetGlobalWorkerNumber( labindex );
            rng( Seed + labindex );
            Set_options_( in_options_ );
        end
    end
end

function Set_options_( in_options_ )
    global options_
    options_ = in_options_;
end
