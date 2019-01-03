function [LIK, LIKK, a, P, rootP] = kalman_filter(Y,start,last,a,P,~,riccati_tol,rescale_prediction_error_covariance,presample,Constant,T,Q,R,H,Z,mm,pp,~,rootF_cond_penalty,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,~,D2P)
% Computes the likelihood of a stationary state space model.

%@info:
%! @deftypefn {Function File} {[@var{LIK},@var{likk},@var{a},@var{P} ] =} DsgeLikelihood (@var{Y}, @var{start}, @var{last}, @var{a}, @var{P}, @var{kalman_tol}, @var{riccati_tol},@var{presample},@var{T},@var{Q},@var{R},@var{H},@var{Z},@var{mm},@var{pp},@var{rr},@var{Zflag},@var{diffuse_periods})
%! @anchor{kalman_filter}
%! @sp 1
%! Computes the likelihood of a stationary state space model, given initial condition for the states (mean and variance).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item Y
%! Matrix (@var{pp}*T) of doubles, data.
%! @item start
%! Integer scalar, first period.
%! @item last
%! Integer scalar, last period (@var{last}-@var{first} has to be inferior to T).
%! @item a
%! Vector (@var{mm}*1) of doubles, initial mean of the state vector.
%! @item P
%! Matrix (@var{mm}*@var{mm}) of doubles, initial covariance matrix of the state vector.
%! @item kalman_tol
%! Double scalar, tolerance parameter (rcond, inversibility of the covariance matrix of the prediction errors).
%! @item riccati_tol
%! Double scalar, tolerance parameter (iteration over the Riccati equation).
%! @item presample
%! Integer scalar, presampling if strictly positive (number of initial iterations to be discarded when evaluating the likelihood).
%! @item T
%! Matrix (@var{mm}*@var{mm}) of doubles, transition matrix of the state equation.
%! @item Q
%! Matrix (@var{rr}*@var{rr}) of doubles, covariance matrix of the structural innovations (noise in the state equation).
%! @item R
%! Matrix (@var{mm}*@var{rr}) of doubles, second matrix of the state equation relating the structural innovations to the state variables.
%! @item H
%! Matrix (@var{pp}*@var{pp}) of doubles, covariance matrix of the measurement errors (if no measurement errors set H as a zero scalar).
%! @item Z
%! Matrix (@var{pp}*@var{mm}) of doubles or vector of integers, matrix relating the states to the observed variables or vector of indices (depending on the value of @var{Zflag}).
%! @item mm
%! Integer scalar, number of state variables.
%! @item pp
%! Integer scalar, number of observed variables.
%! @item rr
%! Integer scalar, number of structural innovations.
%! @item Zflag
%! Integer scalar, equal to 0 if Z is a vector of indices targeting the obseved variables in the state vector, equal to 1 if Z is a @var{pp}*@var{mm} matrix.
%! @item diffuse_periods
%! Integer scalar, number of diffuse filter periods in the initialization step.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item LIK
%! Double scalar, value of (minus) the likelihood.
%! @item likk
%! Column vector of doubles, values of the density of each observation.
%! @item a
%! Vector (@var{mm}*1) of doubles, mean of the state vector at the end of the (sub)sample.
%! @item P
%! Matrix (@var{mm}*@var{mm}) of doubles, covariance of the state vector at the end of the (sub)sample.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{DsgeLikelihood}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{kalman_filter_ss}
%! @end deftypefn
%@eod:

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

if issparse( a ) || issparse( P ) || issparse( T ) || issparse( Q ) || issparse( R ) || issparse( H )
    zerosInternal = @sparse;
else
    zerosInternal = @zeros;
end

% Set defaults.
if nargin<22
    analytic_derivation = 0;
    if nargin<21
        diffuse_periods = 0;
        if nargin<20
            Zflag = 0;
            if nargin<19
                rootF_cond_penalty = 0;
            end
        end
    end
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
rootQ  = robust_root( Q );
rootQQ = R * rootQ;   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
likk = zerosInternal(smpl,1);      % Initialization of the vector gathering the densities.
% LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = true;
asy_hess=0;

if nargout < 5
    rootP = robust_root( P );
else
    rootP = P;
end
rootH = robust_root( H .* ones( size( Y, 1 ) ) );
clear Q H;

if rescale_prediction_error_covariance
    error( 'rescale_prediction_error_covariance is not implemented due to use of square root form.' );
end

if  analytic_derivation == 0
    DLIK=[];
    Hess=[];
    % LIKK=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zerosInternal(k,1);                             % Initialization of the score.
    Da    = zerosInternal(mm,k);                            % Derivative State vector.
    dlikk = zerosInternal(smpl,k);

    % if Zflag==0
    %     C = zerosInternal(pp,mm);
    %     for ii=1:pp, C(ii,Z(ii))=1; end         % SELECTION MATRIX IN MEASUREMENT EQ. (FOR WHEN IT IS NOT CONSTANT)
    % else
    %     C=Z;
    % end
    % dC = zerosInternal(pp,mm,k);   % either selection matrix or schur have zero derivatives
    if analytic_derivation==2
        Hess  = zerosInternal(k,k);                             % Initialization of the Hessian
        D2a    = zerosInternal(mm,k,k);                             % State vector.
        % d2C = zerosInternal(pp,mm,k,k);
    else
        asy_hess=D2T;
        Hess=[];
        D2a=[];
        D2T=[];
        D2Yss=[];
    end
    if asy_hess
        Hess  = zerosInternal(k,k);                             % Initialization of the Hessian
    end
    % LIK={inf,DLIK,Hess};
    % LIKK={likk,dlikk};
end

while notsteady && t<=last
    s = t-start+1;
    if Zflag
        v  = Y(:,t)-Z*a;
        M = qr0( [ rootH.', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * Z.', rootP.' ] );
        % [ G, M ] = qr( [ rootH.', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * Z.', rootP.' ] );
        % M.' * M = M .' * G .' * G * M 
        % = [ rootH.', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP.' * Z.', rootP.' ].' * [ rootH.', zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP.' * Z.', rootP.' ]
        % = [ rootH, Z * rootP; zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP ] * [ rootH.', zerosInternal( size( rootP, 1 ), size( rootH, 2 ) ); rootP.' * Z.', rootP.' ]
        % = [ rootH * rootH.' + Z * rootP * rootP.' * Z.', Z * rootP * rootP.'; rootP * rootP.' * Z.', rootP * rootP.' ]
        % = [ F, F.' * K.'; K * F, P ]
        % = [ M11.', 0; M12.', M22.' ] * [ M11, M12; 0, M22 ] = [ M11.' * M11, M11.' * M12; M12.' * M11, M12.' * M12 + M22.' * M22 ]
        % rootF = M11.'
        % K * rootF * rootF.' = M12.' * M11 = M12.' * rootF.'
        % K * rootF = M12.'
        % P = M12.' * M12 + M22.' * M22 = K * F * K.' + M22.' * M22
        % M22.' * M22 = P - K * F * K.' = P - P * Z.' * iF * F * iF.' * Z * P.' = P - P * Z.' * iF * Z * P
    else
        v  = Y(:,t)-a(Z);
        M = qr0( [ rootH.', zerosInternal( size( rootH, 2 ), size( rootP, 1 ) ); rootP(Z,:).', rootP.' ] );
    end
    rootF = M( 1 : length( v ), 1 : length( v ) ).';
    rootPme = M( ( length( v ) + 1 ) : end, ( length( v ) + 1 ) : end ).';
    K = M( 1 : length( v ), ( length( v ) + 1 ) : end ).' / rootF;

    F_singular = false;
    full_rootF = full( rootF );
    log_dF = sum( log( eig( full_rootF * full_rootF.' ) ) );
    irootFv = rootF \ v;
    likk(s) = log_dF + irootFv.' * irootFv;
    
    if rootF_cond_penalty > 0
        likk(s) = likk(s) + rootF_cond_penalty * log( cond( full_rootF ) ) ^ 4;
    end
    
    M = qr0( [ rootPme.' * T.'; rootQQ.' ] );
    % [ G, M ] = qr( [ rootPme.' * T.'; rootQQ.' ] );
    % M.' * M = M .' * G .' * G * M 
    % = [ rootPme.' * T.'; rootQQ .' ].' * [ rootPme.' * T.'; rootQQ .' ]
    % = [ T * rootPme, rootQQ ] * [ rootPme.' * T.'; rootQQ .' ] 
    % = T * rootPme * rootPme.' * T.' + rootQQ * rootQQ.' 
    % = T * Pme * T.' + QQ
    % = T * ( P - P * Z.' * iF * Z * P ) * T.' + QQ
    rootP = M.';
    
    full_rootP = full( rootP );
    Ptmp = full_rootP * full_rootP.';
    tmp = (a+K*v);
    if analytic_derivation
        iF = inv( full_rootF * full_rootF.' );
        if analytic_derivation==2
            [Da,DP,DLIKt,D2a,D2P, Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2P);
        else
            [Da,DP,DLIKt,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady);
        end
        if t>presample
            DLIK = DLIK + DLIKt;
            if analytic_derivation==2 || asy_hess
                Hess = Hess + Hesst;
            end
        end
        dlikk(s,:)=DLIKt;
    end
    a = Constant + T*tmp;
    P = Ptmp;
    notsteady = max(abs(K(:)-oldK))>riccati_tol;
    oldK = K(:);
    t = t+1;
end

a = full( a );

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Add observation's densities constants and divide by two.
likk(1:s) = .5*(likk(1:s) + pp*log(2*pi));
if analytic_derivation
    DLIK = DLIK/2;
    dlikk = dlikk/2;
    if analytic_derivation==2 || asy_hess
        if asy_hess==0
            Hess = Hess + tril(Hess,-1)';
        end
        Hess = -Hess/2;
    end
end

% Call steady state Kalman filter if needed.
if t <= last
    iF = inv( full_rootF * full_rootF.' );
    if analytic_derivation
        if analytic_derivation==2
            [tmp, tmp2] = kalman_filter_ss(Y, t, last, a, Constant, T, K, iF, log_dF, Z, pp, Zflag, analytic_derivation, Da, DT, DYss, D2a, D2T, D2Yss);
        else
            [tmp, tmp2] = kalman_filter_ss(Y, t, last, a, Constant, T, K, iF, log_dF, Z, pp, Zflag, analytic_derivation, Da, DT, DYss, asy_hess);
        end
        likk(s+1:end) = tmp2{1};
        dlikk(s+1:end,:) = tmp2{2};
        DLIK = DLIK + tmp{2};
        if analytic_derivation==2 || asy_hess
            Hess = Hess + tmp{3};
        end
    else
        [~, likk(s+1:end)] = kalman_filter_ss(Y, t, last, a, Constant, T, K, iF, log_dF, Z, pp, Zflag);
    end
end

% Compute minus the log-likelihood.
if presample>diffuse_periods
    LIK = sum(likk(1+(presample-diffuse_periods):end));
else
    LIK = sum(likk);
end

if analytic_derivation
    if analytic_derivation==2 || asy_hess
        LIK={LIK, DLIK, Hess};
    else
        LIK={LIK, DLIK};
    end
    LIKK={likk, dlikk};
else
    LIKK=likk;
end
