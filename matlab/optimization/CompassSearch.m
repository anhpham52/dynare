function [ x, fx ] = CompassSearch( f, x, lb, ub )

    MaxStep = 0.01;
    PriorStandardDeviation = 0.1;
    PriorStrength = 1;
    TolX = 1e-6;
    TolF = 1e-8;
    RhoScores   = 0.95;
    RhoSteps    = 0.95;
    % RhoGradient = 0.95;
    
    N = length( x );
    D = 2 * N + 2;

    BigFNumber = 1 / TolF;

    InitialStepSizes = min( MaxStep, sqrt( 1 / 12 ) * ( ub - lb ) );
    Directions = [ zeros( N, 2 ), eye( N ), -eye( N ) ];
    StepSizes = [ MaxStep; MaxStep; InitialStepSizes; InitialStepSizes ];
    
    Pool = gcp;
    NumWorkers = Pool.NumWorkers;
    
    NumBlocks = ceil( D / NumWorkers );
    
    MeanScores = zeros( D, 1 );
    MeanScores2 = PriorStandardDeviation * PriorStandardDeviation * ones( D, 1 );
    ScoresSumWeights   = PriorStrength * ones( D, 1 );
    ScoresSumWeights2  = PriorStrength * ones( D, 1 );
    ScoresObservations = PriorStrength * ones( D, 1 );

    MeanGoodSteps  = randn( N, 1 );
    MeanStepsTaken = randn( N, 1 );
    % MeanGradient   = randn( N, 1 );
    % MeanGradientXX = eye( N, N );
    % MeanGradientXY = zeros( N, 1 );

    FMaxChange = 0;

    fx = f( x );

    Iteration = 0;
    
    First = true;
    
    fprintf( '\n' );

    while true
        
        if First && exist( 'CompassSearchState.mat', 'file' )
            load CompassSearchState StepSizes NumWorkers NumBlocks MeanScores MeanScores2 ScoresSumWeights ScoresSumWeights2 ScoresObservations MeanGoodSteps MeanStepsTaken FMaxChange fx x Iteration;
        end
        
        First = false;

        Iteration = Iteration + 1;

        fprintf( '\n\nIteration: %d\n', Iteration );

        VarScores = MeanScores2 - MeanScores .* MeanScores;
        tau = gamrnd( 0.5 * ScoresObservations, 2 ./ max( TolF, ScoresObservations .* VarScores ) );
        ScoreDraw = MeanScores + BigFNumber .* ( StepSizes < TolX ) + randn( D, 1 ) ./ max( TolF, sqrt( tau .* ScoresObservations ) );
        [ ~, Indices ] = sort( ScoreDraw );
        
        Directions( :, 1 ) = MeanGoodSteps ./ norm( MeanGoodSteps );
        Directions( :, 2 ) = MeanStepsTaken ./ norm( MeanStepsTaken );
        % Directions( :, 3 ) = -MeanGradient ./ norm( MeanGradient );

        for Block = 1 : NumBlocks

            BlockIndices = Indices( unique( min( D, ( 1 : NumWorkers ) + ( Block - 1 ) * NumWorkers ) ) );

            if all( StepSizes( BlockIndices ) < TolX )
                break
            end

            BlockLength = length( BlockIndices );

            xNew = zeros( N, BlockLength );
            fNew = zeros( BlockLength, 1 );
            Truncated = false( BlockLength, 1 );

            parfor j = 1 : BlockLength
                i = BlockIndices( j );
                xNew( :, j ) = x + Directions( :, i ) * StepSizes( i ); %#ok<PFBNS>
                if any( xNew( :, j ) < lb ) || any( xNew( :, j ) > ub ) 
                    xNew( :, j ) = max( lb, min( ub, xNew( :, j ) ) );
                    Truncated( j ) = true;
                else
                    Truncated( j ) = false;
                end
                fNew( j ) = f( xNew( :, j ) ); %#ok<PFBNS>
            end

            ScoresSumWeights( BlockIndices )   = RhoScores * ScoresSumWeights( BlockIndices )  + ( 1 - RhoScores );
            ScoresSumWeights2( BlockIndices )  = RhoScores * RhoScores * ScoresSumWeights2( BlockIndices )  + ( 1 - RhoScores ) * ( 1 - RhoScores );
            ScoresObservations( BlockIndices ) = ScoresSumWeights( BlockIndices ) .* ScoresSumWeights( BlockIndices ) ./ ScoresSumWeights2( BlockIndices );
            
            FMaxChange = max( FMaxChange, max( fx - fNew ) );
            
            MeanScores  = max( -FMaxChange, min( FMaxChange, MeanScores ) );
            MeanScores2 = max( -FMaxChange, min( FMaxChange, MeanScores2 ) );
            
            TruncatedFChange = max( -FMaxChange, min( FMaxChange, fNew - fx ) );
            
            MeanScores( BlockIndices )  = RhoScores * MeanScores( BlockIndices )  + ( 1 - RhoScores ) * ( TruncatedFChange      - MeanScores( BlockIndices )  );
            MeanScores2( BlockIndices ) = RhoScores * MeanScores2( BlockIndices ) + ( 1 - RhoScores ) * ( TruncatedFChange .^ 2 - MeanScores2( BlockIndices ) );

            GoodIndices = BlockIndices( fNew < fx );
            BadIndices = BlockIndices( ( fNew >= fx ) | Truncated );

            StepSizes( GoodIndices ) = 1.1 * StepSizes( GoodIndices );
            StepSizes( BadIndices )  = 0.5 * StepSizes( BadIndices );
            
            fprintf( '\nMax step size: %.4g\tMean step size: %.4g', max( StepSizes ), mean( StepSizes ) );
            
            NumNewGoodStepsObservations = sum( fNew < fx );
            if NumNewGoodStepsObservations > 0
                MeanGoodSteps = RhoSteps ^ NumNewGoodStepsObservations * MeanGoodSteps + ( 1 - RhoSteps ^ NumNewGoodStepsObservations ) * ( sum( xNew( :, fNew < fx ) - x, 2 ) - NumNewGoodStepsObservations * MeanGoodSteps );
            end
            
            % DMxNew = xNew - mean( xNew, 2 );
            % DMfNew = fNew - mean( fNew );
            
            % MeanGradientXX = RhoGradient ^ BlockLength * MeanGradientXX + ( 1 - RhoGradient ^ BlockLength ) * ( DMxNew * DMxNew.' - BlockLength * MeanGradientXX );
            % MeanGradientXY = RhoGradient ^ BlockLength * MeanGradientXY + ( 1 - RhoGradient ^ BlockLength ) * ( DMxNew * DMfNew   - BlockLength * MeanGradientXY );
            % MeanGradient   = MeanGradientXX \ MeanGradientXY;

            if ~isempty( GoodIndices )
                [ nfx, BestIndex ] = min( fNew );
                SufficientImprovement = nfx < fx - TolF;
                fx = nfx;
                nx = xNew( :, BestIndex );
                MeanStepsTaken = RhoSteps * MeanStepsTaken + ( 1 - RhoSteps ) * ( nx - x - MeanStepsTaken );
                x = nx;
                fprintf( '\tfx: %.30g', fx );
                if SufficientImprovement
                    break
                end
            end

        end
        
        save CompassSearchState;

        if all( StepSizes < TolX )
            break
        end

    end

end
