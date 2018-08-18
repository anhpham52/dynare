function ParallelInitializationInternal( )
    spmd
        SetGlobalWorkerNumber( labindex );
    end
end
