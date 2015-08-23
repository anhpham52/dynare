function oo_ = convert_dyn_45_to_44(M_, options_, oo_)
%function oo_ = convert_dyn_45_to_44(M_, options_, oo_
% Converts oo_ from 4.5 to 4.4
%
% INPUTS
%    M_          [struct]    dynare model struct
%    options_    [struct]    dynare options struct
%    oo_         [struct]    dynare output struct
%
% OUTPUTS
%    oo_         [struct]    dynare output struct
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2015 Dynare Team
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
% along = with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%% add initial conditions to Bayesian forecasts
if isfield(oo_,'PointForecast')
    var_names=fieldnames(oo_.PointForecast.HPDinf);
    moment_names=fieldnames(oo_.PointForecast);
    for moment_iter=1:length(moment_names)
        for var_iter=1:length(var_names)
            if strcmp(moment_names{moment_iter},'deciles')
                oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(:,end)*ones(M_.maximum_endo_lag,1)  oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})];
                oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(:,end)*ones(M_.maximum_endo_lag,1) oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})];                
            else
                oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(end)*ones(M_.maximum_endo_lag,1); oo_.MeanForecast.(moment_names{moment_iter}).(var_names{var_iter})];
                oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})=...
                    [oo_.SmoothedVariables.(moment_names{moment_iter}).(var_names{var_iter})(end)*ones(M_.maximum_endo_lag,1); oo_.PointForecast.(moment_names{moment_iter}).(var_names{var_iter})];
            end
        end
    end
end

%% change HPD-fields back to row vectors
if isfield(oo_.PointForecast,'HPDinf')
    names=fieldnames(oo_.PointForecast.HPDinf);
    for ii=1:length(names)
        oo_.PointForecast.HPDinf.(names{ii})=oo_.PointForecast.HPDinf.(names{ii})';
        oo_.PointForecast.HPDsup.(names{ii})=oo_.PointForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_.MeanForecast,'HPDinf')
    names=fieldnames(oo_.MeanForecast.HPDinf);
    for ii=1:length(names)
        oo_.MeanForecast.HPDinf.(names{ii})=oo_.MeanForecast.HPDinf.(names{ii})';
        oo_.MeanForecast.HPDsup.(names{ii})=oo_.MeanForecast.HPDsup.(names{ii})';
    end
end

if isfield(oo_.UpdatedVariables,'HPDinf')
    names=fieldnames(oo_.UpdatedVariables.HPDinf);
    for ii=1:length(names)
        oo_.UpdatedVariables.HPDinf.(names{ii})=oo_.UpdatedVariables.HPDinf.(names{ii})';
        oo_.UpdatedVariables.HPDsup.(names{ii})=oo_.UpdatedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.SmoothedVariables,'HPDinf')
    names=fieldnames(oo_.SmoothedVariables.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedVariables.HPDinf.(names{ii})=oo_.SmoothedVariables.HPDinf.(names{ii})';
        oo_.SmoothedVariables.HPDsup.(names{ii})=oo_.SmoothedVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.FilteredVariables,'HPDinf')
    names=fieldnames(oo_.FilteredVariables.HPDinf);
    for ii=1:length(names)
        oo_.FilteredVariables.HPDinf.(names{ii})=oo_.FilteredVariables.HPDinf.(names{ii})';
        oo_.FilteredVariables.HPDsup.(names{ii})=oo_.FilteredVariables.HPDsup.(names{ii})';
    end
end

if isfield(oo_.SmoothedShocks,'HPDinf')
    names=fieldnames(oo_.SmoothedShocks.HPDinf);
    for ii=1:length(names)
        oo_.SmoothedShocks.HPDinf.(names{ii})=oo_.SmoothedShocks.HPDinf.(names{ii})';
        oo_.SmoothedShocks.HPDsup.(names{ii})=oo_.SmoothedShocks.HPDsup.(names{ii})';
    end
end

%% padd classical filtered variables with redundant zeros
if isfield(oo_,'FilteredVariables')
    names=fieldnames(oo_.FilteredVariables);
    for ii=1:length(names)
        %make sure Bayesian fields are not affect
        if ~strcmp(names{ii},'Mean') && ~strcmp(names{ii},'Median') && ~strcmp(names{ii},'deciles') ...
                && ~strcmp(names{ii},'Var') && ~strcmp(names{ii},'HPDinf') && ~strcmp(names{ii},'HPDsup') 
            oo_.FilteredVariables.(names{ii})=[0; oo_.FilteredVariables.(names{ii}); zeros(options_.nk-1,1)];
        end
    end
end

%% set old field posterior_std and remove new field posterior_std_at_mode
if isfield(oo_,'posterior_std_at_mode')
    oo_.posterior_std=oo_.posterior_std_at_mode;
    oo_=rmfield(oo_,'posterior_std_at_mode');
end
