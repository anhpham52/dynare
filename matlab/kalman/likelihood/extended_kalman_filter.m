function  [LIK, lik, a, P, rootP] = extended_kalman_filter(data_index,~,~,Y,start,last,a,P,~,~,rescale_prediction_error_covariance,presample,StateSelect,Constant,Jacobian0,Hessian,Q,H,Z,~,pp,~,Zflag,diffuse_periods)
% Computes the likelihood of a state space model in the case with missing observations.
%
% INPUTS
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]    pp*smpl matrix of data.
%    start                        [integer]   scalar, index of the first observation.
%    last                         [integer]   scalar, index of the last observation.
%    a                            [double]    pp*1 vector, initial level of the state vector.
%    P                            [double]    pp*pp matrix, covariance matrix of the initial state vector.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).
%    presample                    [integer]   scalar, presampling if strictly positive.
%    T                            [double]    mm*mm transition matrix of the state equation.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    H                            [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors.
%    Z                            [integer]   pp*1 vector of indices for the observed variables.
%    mm                           [integer]   scalar, dimension of the state vector.
%    pp                           [integer]   scalar, number of observed variables.
%    rr                           [integer]   scalar, number of structural innovations.
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%    a          [double]    mm*1 vector, estimated level of the states.
%    P          [double]    mm*mm matrix, covariance matrix of the states.
%
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if issparse( a ) || issparse( P ) || issparse( Constant ) || issparse( Q ) || issparse( Jacobian0 ) || issparse( H ) || issparse( Hessian )
    zerosInternal = @sparse;
else
    zerosInternal = @zeros;
end

% Set defaults
if nargin<20
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<21
    diffuse_periods = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

if isequal(H,0)
    H = zerosInternal(pp,pp);
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
%dF   = 1;
% rootQ  = robust_root( Q );
% rootQQ = R * rootQ;   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zerosInternal(smpl,1);      % Initialization of the vector gathering the densities.
%LIK  = Inf;                % Default value of the log likelihood.
%oldK = Inf;
notsteady   = 1;
F_singular  = true;
s = 0;

nState = length( StateSelect );
StateBlock = 1 : nState;
na = length( a );

if nargout < 5
    rootP = robust_root( P );
else
    rootP = P;
end

rootH = robust_root( H .* ones( size( Y, 1 ) ) );
clear H;

if rescale_prediction_error_covariance
    error( 'rescale_prediction_error_covariance is not implemented due to use of square root form.' );
end

while notsteady && t<=last
    s  = t-start+1;
    d_index = data_index{t};

    aHat = a( StateSelect );
    
    P = rootP * rootP.';

    Ptilde = blkdiag( P( StateSelect, StateSelect ), Q );
    
    a = Constant + Jacobian0( :, StateBlock ) * aHat;
    
    Jacobian = Jacobian0;
    
    PtildeHessian = cell( na, 1 );
    
    for j = 1 : na
        PtildeHessian{ j } = Ptilde * Hessian{ j };
        aHatHessian = aHat .' * Hessian{ j }( StateBlock, : );
        
        a( j ) = a( j ) + 0.5 * ( aHatHessian( :, StateBlock ) * aHat ) + 0.5 * trace( PtildeHessian{ j } );
        Jacobian( j, : ) = Jacobian( j, : ) + aHatHessian;
    end
    
    P = Jacobian * Ptilde * Jacobian.';
    
    for j = 1 : na
        for i = j : na
            P( i, j ) = P( i, j ) + 0.5 * sum( PtildeHessian{ i }(:) .* PtildeHessian{ j }(:) );
            P( j, i ) = P( i, j );
        end
    end
    
    rootP = robust_root( P );
       
    if ~isempty(d_index)
        % Compute the prediction error and its variance
        if Zflag
            z = Z(d_index,:);
            v = Y(d_index,t)-z*a;
            M = qr0( [ rootH( d_index, : ).', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * z.', rootP.' ] );
            % [ G, M ] = qr( [ rootH( d_index, : ).', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * z.', rootP.' ] );
            % M.' * M = M .' * G .' * G * M 
            % = [ rootH( d_index, : ).', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * z.', rootP.' ].' * [ rootH( d_index, : ).', zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP.' * z.', rootP.' ]
            % = [ rootH( d_index, : ), z * rootP; zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP ] * [ rootH( d_index, : ).', zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP.' * z.', rootP.' ]
            % = [ rootH( d_index, : ) * rootH( d_index, : ).' + z * rootP * rootP.' * z.', z * rootP * rootP.'; rootP * rootP.' * z.', rootP * rootP.' ]
            % = [ F, F.' * K.'; K * F, P ]
            % = [ M11.', 0; M12.', M22.' ] * [ M11, M12; 0, M22 ] = [ M11.' * M11, M11.' * M12; M12.' * M11, M12.' * M12 + M22.' * M22 ]
            % rootF = M11.'
            % K * rootF * rootF.' = M12.' * M11 = M12.' * rootF.'
            % K * rootF = M12.'
            % P = M12.' * M12 + M22.' * M22 = K * F * K.' + M22.' * M22
            % M22.' * M22 = P - K * F * K.' = P - P * z.' * iF * F * iF.' * z * P.' = P - P * z.' * iF * z * P
        else
            z = Z(d_index);
            v = Y(d_index,t) - a(z);
            M = qr0( [ rootH( d_index, : ).', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP(z,:).', rootP.' ] );
        end
        rootF = M( 1 : length( d_index ), 1 : length( d_index ) ).';
        rootP = M( ( length( d_index ) + 1 ) : end, ( length( d_index ) + 1 ) : end ).';
        K = M( 1 : length( d_index ), ( length( d_index ) + 1 ) : end ).' / rootF;

        F_singular = false;
        full_rootF = full( rootF );
        log_dF = sum( log( eig( full_rootF * full_rootF.' ) ) );
        irootFv = rootF \ v;
        lik(s) = log_dF + irootFv.' * irootFv + length(d_index)*log(2*pi);
        
        a = a + K * v;
    end
    t = t+1;
end

full_rootP = full( rootP );
P = full_rootP * full_rootP.';

a = full( a );

if F_singular && smpl > 1
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Divide by two.
lik(1:s) = .5*lik(1:s);

% Compute minus the log-likelihood.
if presample>=diffuse_periods
    LIK = sum(lik(1+presample-diffuse_periods:end));
else
    LIK = sum(lik);
end
