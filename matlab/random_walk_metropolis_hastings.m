function random_walk_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% function random_walk_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% Random Walk Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o ProposalFun  [char]   string specifying the name of the proposal
%                           density
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o dataset_              data structure
%   o dataset_info          dataset info structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%
% ALGORITHM 
%   Random-Walk Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code suitable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle)
% has been removed from this function and been placed in the random_walk_metropolis_hastings_core.m funtion.
% 
% The DYNARE parallel packages comprise a i) set of pairs of Matlab functions that can be executed in
% parallel and called name_function.m and name_function_core.m and ii) a second set of functions used
% to manage the parallel computations.
%
% This function was the first function to be parallelized. Later, other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management functions.

% Copyright (C) 2006-2015 Dynare Team
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


% In Metropolis, we set penalty to Inf so as to reject all parameter sets triggering an error during target density computation
global objective_function_penalty_base
objective_function_penalty_base = Inf;

% Initialization of the random walk metropolis-hastings chains.
[ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
    metropolis_hastings_initialization(TargetFun, xparam1, vv, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_);

InitSizeArray = min([repmat(MAX_nruns,nblck,1) fline+nruns-1],[],2);

% Load last mh history file
load_last_mh_history_file(MetropolisFolder, ModelName);

% Only for test parallel results!!!

% To check the equivalence between parallel and serial computation!
% First run in serial mode, and then comment the follow line.
%   save('recordSerial.mat','-struct', 'record');

% For parallel runs after serial runs with the abobe line active.
%   TempRecord=load('recordSerial.mat');
%   record.Seeds=TempRecord.Seeds;



% Snapshot of the current state of computing. It necessary for the parallel
% execution (i.e. to execute in a corretct way a portion of code remotely or
% on many cores). The mandatory variables for local/remote parallel
% computing are stored in the localVars struct.

if options_.TaRB.use_TaRB
    options_.silent_optimizer=1; %locally set optimizer to silent mode
end

localVars =   struct('TargetFun', TargetFun, ...
                     'ProposalFun', ProposalFun, ...
                     'xparam1', xparam1, ...
                     'vv', vv, ...
                     'mh_bounds', mh_bounds, ...
                     'ix2', ix2, ...
                     'ilogpo2', ilogpo2, ...
                     'ModelName', ModelName, ...
                     'fline', fline, ...
                     'npar', npar, ...
                     'nruns', nruns, ...
                     'NewFile', NewFile, ...
                     'MAX_nruns', MAX_nruns, ...
                     'd', d, ...
                     'InitSizeArray',InitSizeArray, ...
                     'record', record, ...
                     'dataset_', dataset_, ...
                     'dataset_info', dataset_info, ...
                     'options_', options_, ...
                     'M_',M_, ...
                     'bayestopt_', bayestopt_, ...
                     'estim_params_', estim_params_, ...
                     'oo_', oo_,...
                     'varargin',[]);


% User doesn't want to use parallel computing, or wants to compute a
% single chain compute Random walk Metropolis-Hastings algorithm sequentially.

if isnumeric(options_.parallel) || (nblck-fblck)==0,
    if options_.TaRB.use_TaRB
        fout = TaRB_metropolis_hastings_core(localVars, fblck, nblck, 0);
    else
        fout = random_walk_metropolis_hastings_core(localVars, fblck, nblck, 0);
    end
    record = fout.record;
    % Parallel in Local or remote machine.   
else 
    % Global variables for parallel routines.
    globalVars = struct();
    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[ModelName '_static.m']};
    NamFileInput(2,:) = {'',[ModelName '_dynamic.m']};
    if options_.steadystate_flag,
        NamFileInput(length(NamFileInput)+1,:)={'',[ModelName '_steadystate.m']};
    end
    if (options_.load_mh_file~=0)  && any(fline>1) ,
        NamFileInput(length(NamFileInput)+1,:)={[M_.dname '/metropolis/'],[ModelName '_mh' int2str(NewFile(1)) '_blck*.mat']};
    end
    if exist([ModelName '_optimal_mh_scale_parameter.mat'])
        NamFileInput(length(NamFileInput)+1,:)={'',[ModelName '_optimal_mh_scale_parameter.mat']};
    end
    % from where to get back results
    %     NamFileOutput(1,:) = {[M_.dname,'/metropolis/'],'*.*'};
    if options_.TaRB.use_TaRB
        [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, fblck, nblck,NamFileInput,'TaRB_metropolis_hastings_core', localVars, globalVars, options_.parallel_info);        
    else    
        [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, fblck, nblck,NamFileInput,'random_walk_metropolis_hastings_core', localVars, globalVars, options_.parallel_info);
    end
    for j=1:totCPU,
        offset = sum(nBlockPerCPU(1:j-1))+fblck-1;
        record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)));
        record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:)=fout(j).record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:);
        record.AcceptanceRatio(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.AcceptanceRatio(offset+1:sum(nBlockPerCPU(1:j)));
        record.LastSeeds(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.LastSeeds(offset+1:sum(nBlockPerCPU(1:j)));
    end

end

irun = fout(1).irun;
NewFile = fout(1).NewFile;

record.MCMCConcludedSuccessfully = 1; %set indicator for successful run

update_last_mh_history_file(MetropolisFolder, ModelName, record);

% Provide diagnostic output
skipline()
disp(['Estimation::mcmc: Number of mh files: ' int2str(NewFile(1)) ' per block.'])
disp(['Estimation::mcmc: Total number of generated files: ' int2str(NewFile(1)*nblck) '.'])
disp(['Estimation::mcmc: Total number of iterations: ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp(['Estimation::mcmc: Current acceptance ratio per chain: '])
for i=1:nblck
    if i<10
        disp(['                                                       Chain  ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    else
        disp(['                                                       Chain ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    end
end