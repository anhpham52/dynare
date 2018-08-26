function ParallelInitializationInternal( varargin )
    if nargin == 0
        spmd
            SetGlobalWorkerNumber( labindex );
        end
    else
        in_options_ = varargin{ 1 };
        spmd
            SetGlobalWorkerNumber( labindex );
            Set_options_( in_options_ );
        end
    end
end

function Set_options_( in_options_ )
    global options_
    options_ = in_options_;
end
