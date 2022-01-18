%% mnb_PSOFit
% Run the particle swarm optimization routine.
%
% Usage: 
%
% * [gp_best,gf_best,ExitReason] = mnb_PSOFit(ObjFunName,LB,UB)
% * [gp_best,gf_best,ExitReason] = mnb_PSOFit(ObjFunName,LB,UB,varargin)
% * [gp_best,gf_best,ExitReason] = mnb_PSOFit(ObjFunName,LB,UB,PSOoptsStruct)
%

%% Required Input Arguments (in order)
%
% * ObjFuncName : char : name of the user-supplied objective function
%
% * LB : numeric vector : lower bounds for estimated parameter values. 
%
% * UB : numeric vector : upper bounds for estimated parameter values. 
%
%% Optional Input Arguments (passed via Param-Value pairs)
%
% * ObjFuncDir : char : path to the directory containing the user-supplied
%       objective function. Default = pwd
%
% * ObjFuncData : any variable type : The argument is simply
%       passed-through from this fitting function to the objective
%       function. It can be used, for example, to store the experimental
%       data, a SimBiology model, etc.
%
% * InitialGuess : numeric vector of the same length as LB and UB
%       arguments : These initial parameter values will be used to
%       initialize the position of one of the starting particles. An error
%       will be thrown if any of these values are outside of the range
%       defined by LB and UB.
%
% * LogTransform : logical vector of the same length as LB and UB : this
%       vector defines which parameters will be log-transformed (true) or
%       not (false). Log transformation is a useful technique when a
%       parameter's boundaries span 2 or more orders of magnitude. The
%       values will be linearized when passed to the objective function.
%       Default = [] (do not log-transform any parameters)
%
% * Postions : numeric array of size NumParameters x NumSwp : This input
%       argument contains the positions matrix from a prior PSO fitting run
%       and is intended to allow the user to restart a new fit from that
%       prior state. It must be accompanied by the Velocities input
%       argument (below). Default = []
%
% * Velocities : numeric array of size NumParameters x NumSwp : This input
%       argument contains the velocities matrix from a prior PSO fitting
%       run and is intended to allow the user to restart a new fit from
%       that prior state. It must be accompanied by the Positions input
%       argument (above). Default = []
%
% * Verbose : logical : display status updates. Default = true
%
% * FN_OptState : char : MAT file to which the final fitness value and
%       parameter values are saved. This file will also contain the full
%       positions and velocities of every particle in the swarm as well as
%       all PSO status parameters (loop variable, count, swarm size, etc). 
%       Default = 'PSO Fit - Final Results.mat'
%
% * FN_IntState : char : MAT file to which the interim fitting results
%       are saved. This file will contain the same information as in the
%       "FN_OptState" file; the difference being that this information is
%       from an interim fit state.
%       Default = 'PSO Fit - interim results.mat'
%
% * FN_IntState_restart : char : MAT file name with the interim fit results
%       from a previous PSO run. It is in the same format as "FN_IntState".
%       This option enables a prior fit routine to be restarted. 
%       Default = '' (do not load past results)
%
% * MaxIter : positive integer >= 1 : maximum number of iterations that can
%       be performed in the PSO loop. Default = 2500
%
% * FitCount : numeric scalar > 0 : The PSO loop can exit if the "count"
%       variable (running tally of non-improving iterations) reaches this
%       limit. Default = 100.
%
% * NumSwp : positive integer : Number of particles in the swarm. 
%       Default = 20
%
% * wmin : numeric scalar > 0 : Minimum intertia factor. Default = 0.1
%
% * wmax : numeric scalar > 0 : Maximum intertia factor. Default = 1.1
%
% * wfactor_up : numeric scalar > 0 : When the 
%
% * wfactor_dn : numeric scalar > 0 : 
%
% * wcount_min : integer >=1 : Counter value (for non-improving iterations)
%       where the value of the adaptive inertia factor begins to increase.
%       Default = 4
%
% * wcount_max : integer > wcount_min : Counter value (for non-improving
%       iterations) where the value of the adaptive inertia factor begins 
%       to decrease. Default = 7
%
% * c1 : numeric scalar > 0 : Self-recognition and social component
%       coefficient 1. Default = 1.49
% 
% * c2 : numeric scalar > 0 : Self-recognition and social component
%       coefficient 2. Default = 1.49
%
% * lmin : integer > 1 and < NumSwp : Randomly assigned neighborhoods
%       Default = 5
%
% * InitRandStream : logical : Initialize random generator
%
% * UseParallel : logical : Use to run parallelized jobs on clusters
%
% * ParticleElapsedTime : logical : Compute time elapsed by each particle
% at any given PSO iteration
%
%% Function Output Arguments
%
% * gp_best : numeric vector of the best parameter values
%
% * gf_best : numeric scalar: fitness of the best parameter values
%
% * ExitReason : char : Reason for PSO exit. One of the following messages
% will be returned... (1) "MaxIter Reached" (2) "FitCount Reached"
%
%% Dependencies
%
% * mnb_PSOOptions -- required to get default input argument values.
%
%% Limitations
%
% * TODO: Can we get a speed boost anywhere? Vectorizing more loops?  
%
% * TODO: add more informative error messages, especially to the PSO
% restart code
%
%% Version Information
% * Last changed date: $Date: 2015-03-04 10:49:38 -0500 (Wed, 04 Mar 2015) $
% * Last changed by: $Author: jkearns $
% * Last changed revision: $Revision: 1560 $
%
% Copyright 2012-2015 Merrimack Pharmaceuticals, Inc.

function [gp_best,gf_best,ExitReason] = mnb_PSOFit(ObjFuncName,LB,UB,varargin)

%-- Save the clock time for later reporting of total elapsed time
    StartClock = clock;

%-- Validate and parse the input

    % Get Default PSO options
    DefPSO = mnb_PSOOptions();

    % Configure the Parser
    p = inputParser();
    p.FunctionName  = 'MNB_PSOFIT';
    p.CaseSensitive = false;
    p.StructExpand  = true;
    p.KeepUnmatched = true; %ignore unknown input arguments

    % Validate Required Input Arguments (must be passed in this order)
    %p.addRequired('m1',@(x) isa(x,'SimBiology.Model'));
    p.addRequired('ObjFuncName',@(x) ischar(x) && ~isempty(x)); 
    p.addRequired('LB',@(x) isnumeric(x) && ~isempty(x) && (isrow(x) || iscolumn(x)) && (numel(x) == numel(UB)));
    p.addRequired('UB',@(x) isnumeric(x) && ~isempty(x) && (isrow(x) || iscolumn(x)) && (numel(x) == numel(LB)));

    % Validate Optional Input Arguments (use Param-Value pairs)
    p.addParameter('ObjFuncDir',pwd,@(x) isempty(x) || (ischar(x) && exist(x, 'dir')));
    p.addParameter('ObjFuncData',@(x) ~isempty(x)); % any variable type

    p.addParameter('InitialGuess',DefPSO.InitialGuess,@(x) isempty(x) || (isnumeric(x) && isvector(x) && (numel(x) == numel(LB)))); % empty or numeric vector of same length as LB
    p.addParameter('LogTransform',[],@(x) isempty(x) || (islogical(x) && isvector(x) && (numel(x) == numel(LB)))); % empty or logical vector of same length as LB

    p.addParameter('NumSwp', DefPSO.NumSwp, @(x) isnumeric(x) && isscalar(x) && (x>=1) && (mod(x,1) == 0)); % positive integer
    p.addParameter('lmin', DefPSO.lmin, @(x) isnumeric(x) && isscalar(x) && (x>=1) && (mod(x,1) == 0)); % positive integer
    p.addParameter('c1', DefPSO.c1, @(x) isnumeric(x) && isscalar(x) && (x>0));        % positive real number
    p.addParameter('c2', DefPSO.c2, @(x) isnumeric(x) && isscalar(x) && (x>0));        % positive real number
    
    p.addParameter('wmin', DefPSO.wmin, @(x) isnumeric(x) && isscalar(x) && (x>0));    % positive real number
    p.addParameter('wmax', DefPSO.wmax, @(x) isnumeric(x) && isscalar(x) && (x>0));    % positive real number
    p.addParameter('wcount_min', DefPSO.wcount_min, @(x) isnumeric(x) && isscalar(x) && (x>=1) && (mod(x,1) == 0)); % positive integer
    p.addParameter('wcount_max', DefPSO.wcount_max, @(x) isnumeric(x) && isscalar(x) && (x>=1) && (mod(x,1) == 0)); % positive integer
    p.addParameter('wfactor_up', DefPSO.wfactor_up, @(x) isnumeric(x) && isscalar(x) && (x>0)); % positive real number
    p.addParameter('wfactor_dn', DefPSO.wfactor_dn, @(x) isnumeric(x) && isscalar(x) && (x>0)); % positive real number

    p.addParameter('FitCount',100,@(x) isnumeric(x) && isscalar(x) && (x>=1));
    p.addParameter('MaxIter', DefPSO.MaxIter, @(x) isnumeric(x) && isscalar(x) && (x>=1) && (mod(x,1) == 0)); % positive integer

    p.addParameter('Verbose',false,@islogical);
    p.addParameter('FN_OptState',DefPSO.FN_OptState,@(x) isempty(x) || ischar(x));
    p.addParameter('FN_IntState',DefPSO.FN_IntState,@(x) isempty(x) || ischar(x));
    p.addParameter('FN_IntState_restart',DefPSO.FN_IntState_restart,@(x) isempty(x) || ischar(x));
    p.addParameter('ExtendMaxIter',false,@islogical);

    p.addParameter('InitRandStream',DefPSO.InitRandStream,@islogical); % scalar logical
    p.addParameter('UseParallel',DefPSO.UseParallel,@islogical); % scalar logical
    p.addParameter('ParticleElapsedTime',DefPSO.ParticleElapsedTime,@islogical); % scalar logical

    % Parse the input arguments and push to the struct 'r'
    p.parse(ObjFuncName,LB,UB,varargin{:});  
    r = p.Results();

    % Post-processing of parsed input arguments (mostly error checking)

        % Convert LB, UB, InitialGuess, and LogTransform vectors to column (if row)
        if isrow(r.LB); r.LB = r.LB'; end
        if isrow(r.UB); r.UB = r.UB'; end

        if ~isempty(r.InitialGuess) && isrow(r.InitialGuess) 
            r.InitialGuess = r.InitialGuess'; 
        end  

        if ~isempty(r.LogTransform) && isrow(r.LogTransform) 
            r.LogTransform = r.LogTransform'; 
        end  

        % Confirm that UB values are >= to the corresponding LB values
        if ~all(r.UB-r.LB >= 0)
            error([p.FunctionName ':BoundValues'],...
               ['All Upper Bound values must be >= to the '...
               'corresponding Lower Bound values.']);
        end

        % Confirm that all InitialGuess values are between the bounds
        if ~isempty(r.InitialGuess) && (~all(r.InitialGuess <= r.UB) || ~all(r.InitialGuess >= r.LB))
            error([p.FunctionName ':InitialGuessOutOfBounds'],...
               'All InitialGuess values must be within the LB and UB values.');
        end

        % Log transform UB, LB, and InitialGuess values
        if ~isempty(r.LogTransform)

            % Throw an error if bound values are <= 0
            if any(r.LB(r.LogTransform) <= 0)
                error([p.FunctionName ':LogTransformLB'],...
                'Cannot log transform LB. One or more values are <= 0');
            end

            if any(r.UB(r.LogTransform) <= 0)
                error([p.FunctionName ':LogTransformUB'],...
                'Cannot log transform UB. One or more values are <= 0');
            end

            r.LB(r.LogTransform) = log10(r.LB(r.LogTransform));
            r.UB(r.LogTransform) = log10(r.UB(r.LogTransform));

            if ~isempty(r.InitialGuess)

                if any(r.InitialGuess(r.LogTransform) <= 0)
                    error([p.FunctionName ':LogTransformUB'],...
                    'Cannot log transform InitialGuess. One or more values are <= 0');
                end
                r.InitialGuess(r.LogTransform) = log10(r.InitialGuess(r.LogTransform));
            end

            UseLogTransform = true;
            LogTransIdx    = r.LogTransform; % Needed for saving results to a MAT file (below)

        else
            UseLogTransform = false;
            LogTransIdx    = false(length(r.LB),1);
        end

        % Confirm that lmin is <= NumSwp-1
        %   This says that the neighborhood size for a particle cannot
        %   be greater than the number of the other particles
        if r.lmin > r.NumSwp
            error([p.FunctionName ':lminOutOfBounds'],...
               'The value of lmin must be <= NumSwp');
        end

        % Confirm that count_max value is > than count_min value
        if r.wcount_max <= r.wcount_min
            error(['maximum value of counter for adjusting inertia ' ...
                'factor is bigger than minimum value of counter.' ...
                '(wcount_max must be > wcount_min)']);
        end

        % Convert the objective function name to a function handle
        if ~isempty(r.ObjFuncDir)
            p = path; % Save the MATLAB path
            path(r.ObjFuncDir,path);
            ObjFuncHandle = str2func(r.ObjFuncName);
            path(p)
        else
            ObjFuncHandle = str2func(r.ObjFuncName);
        end

        % Initialize random generator
        if r.InitRandStream  % set to false if randstream is initialized externally
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
        end

%-- Initialize the swarm
%       Option 1: Initialize a new fit
%       Option 2: Restart from a prior PSO routine

    NumPar      = numel(r.LB);

    if isempty(r.FN_IntState_restart) % Option 1: Initialize a new fit

        % Initialize swarm particles and model parameters dimension
            lp_best         = zeros(NumPar,r.NumSwp); % best positions for each particle
            gp_best         = zeros(NumPar,1);        % global best positions 

        % Initialize counters
            StartIter = 1; count = 0;   l = r.lmin; 

        % Intialize the inertia factor, local best vector and global best
            w = r.wmax ;   lf_best = inf*ones(r.NumSwp,1);   gf_best = inf;  

        % Intialize cutoff criteria
            gf_bestHis      = nan(r.MaxIter,1); 

        % Create a position matrix of size (NumPar,r.NumSwp) with randomly
        % assigned positions between the LB and UB bound values.
        %  This is of the form: LB + rand_matrix*(UB-LB)

            pos = bsxfun(@plus,r.LB,bsxfun(@times,r.UB-r.LB,rand(NumPar,r.NumSwp)));

            % If the user supplied an initial guess, assign it to the first
            % positions column (param values for the first particle)
            if ~isempty(r.InitialGuess)
                pos(:,1) = r.InitialGuess;
            end

            % Set the starting velocities to zero
            vel = zeros(NumPar,r.NumSwp);
            
        % Initialize matrix for particle elapsed time
        if r.ParticleElapsedTime
            SwpElapsedTime = nan(r.NumSwp,r.MaxIter);
        else
            SwpElapsedTime = [];
        end         

    else % Option 2: restart a fit from a previous state

        % Read in the file
            LoadVars = load(r.FN_IntState_restart);

        % Confirm that all variables are present
            StateVars = {'pos','vel','k','count','l','w',...
                     'lf_best','gf_best','lp_best','gp_best','gf_bestHis' 'SwpElapsedTime'};

            CheckVars = isfield(LoadVars,StateVars);

            if ~all(CheckVars)
                error([p.FunctionName ':IntState_restart1'],...
                    ['Not all fit variables are present in '...
                    'the supplied MAT file.']);
            end

        % Confirm that the size of the positions and velocities
        % matrices are consistent with the options that the user
        % has set (#params, #particles).
            [posCol,posRow] = size(LoadVars.pos);
            if ~isequal(posCol,NumPar) || ~isequal(posRow,r.NumSwp)
                disp('hi1');
                error([p.FunctionName ':IntState_restart2'],...
                      ['The size of the positions matrix must be ',...
                      ' NumParams (%d) x NumSwp (%d).'],NumPar,r.NumSwp);
            end

            [velCol,velRow] = size(LoadVars.vel);
            if ~isequal(velCol,NumPar) || ~isequal(velRow,r.NumSwp)
                disp('hi2');
                error([p.FunctionName ':IntState_restart3'],...
                      ['The size of the velocities matrix must be ',...
                      ' NumParams (%d) x NumSwp (%d).'],NumPar,r.NumSwp);
            end

        % Confirm that all positions are within the bounds (LB,UB)
            for np = 1:NumPar
                if ~all(LoadVars.pos(np,:) >= r.LB(np))
                    error([p.FunctionName ':IntState_restart4'],...
                        ['One or more values in the supplied Positions '...
                        'matrix are below lower bound #%d'],np);
                end

                if ~all(LoadVars.pos(np,:) <= r.UB(np))
                    error([p.FunctionName ':IntState_restart5'],...
                        ['One or more values in the supplied Positions '...
                        'matrix are above upper bound #%d'],np); 
                end
            end

        % Initialize the swarm particle positions and velocities
            pos = LoadVars.pos;
            vel = LoadVars.vel;

        % Initialize swarm particles and model parameters dimension
            lp_best         = LoadVars.lp_best; % best positions for each particle
            gp_best         = LoadVars.gp_best; % global best positions 

        % Initialize counters
            StartIter = LoadVars.k; count = LoadVars.count;   l = LoadVars.l; 

        % Intialize the inertia factor, local best vector and global best
            w = LoadVars.w ;   lf_best = LoadVars.lf_best;   gf_best = LoadVars.gf_best;  

        % Intialize cutoff criteria
            gf_bestHis      = LoadVars.gf_bestHis;
            % Expand the size of gf_bestHis if restarting/continuing a simulation that hit MaxIterNum
            if r.ExtendMaxIter
                if length(gf_bestHis)< r.MaxIter
                    gf_bestHis = [gf_bestHis; NaN(r.MaxIter-length(gf_bestHis),1)];
                end
            end
        % Initialize matrix for particle elapsed time
            % Expand the size of SwpElapsedTime if restarting/continuing a simulation that hit MaxIterNum
            SwpElapsedTime = LoadVars.SwpElapsedTime;
            if r.ExtendMaxIter
                if size(SwpElapsedTime,2)< r.MaxIter
                    SwpElapsedTime = [SwpElapsedTime zeros(size(SwpElapsedTime,1),r.MaxIter-size(SwpElapsedTime,2))];
                end
            end
            % Confirm that size is correct
            if ~isequal(size(SwpElapsedTime),[r.NumSwp r.MaxIter])
                error([p.FunctionName ':IntState_restart6'],...
                        ['The size of "SwpElapsedTime" must equal '...
                        '#Swp(%d) x #MaxIter(%d)'],r.NumSwp,r.MaxIter);
            end

        % For good measure, clear the LoadVars variable from the workspace
            clear LoadVars;

    end

    % Initialize other options that are in common for restart or fresh fits

        % vector of fitness values returned from objective function
        local_fitness   = zeros(r.NumSwp,1);

        % Grab ObjFuncData out of the parser struct to better enable the parfor
        ObjFuncData     = r.ObjFuncData; 

        % Build a matrix to store swarm particle numbers (used in the selection
        % of a local neighborhood within the MaxIter loop. The vector is sized
        % r.NumSwp-1 by r.NumSwp. Each column contains the index of each
        % particle in the swarm except for that particle.
        SwmNeigh_noSelf = zeros(r.NumSwp-1,r.NumSwp);
        swarm_vec       = (1:r.NumSwp)';

        for np = 1:r.NumSwp
            SwmNeigh_noSelf(:,np) = swarm_vec(swarm_vec ~= np);
        end

        % Temporary variable for final print on screen
        BreakFlag = false;

%-- Swarm Optimization Loop

    for k = StartIter:r.MaxIter

        % Evaluate fitness of each swarm particle

            % Linearize any log transformed parameter values before passing
            % to the objective function
            if UseLogTransform
                LinPos = pos;
                LinPos(LogTransIdx,:) = 10.^(LinPos(LogTransIdx,:));

                if r.UseParallel
                    if r.ParticleElapsedTime
                        parfor ns = 1:r.NumSwp
                            ParticleStart = clock;
                            local_fitness(ns) = ObjFuncHandle(LinPos(:,ns),ObjFuncData); 
                            SwpElapsedTime(ns,k) = etime(clock,ParticleStart);
                        end
                    else
                        parfor ns = 1:r.NumSwp
                            local_fitness(ns) = ObjFuncHandle(LinPos(:,ns),ObjFuncData); 
                        end
                    end

                else
                    for ns = 1:r.NumSwp
                        if r.ParticleElapsedTime
                            ParticleStart = clock;
                            local_fitness(ns) = ObjFuncHandle(LinPos(:,ns),ObjFuncData); 
                            SwpElapsedTime(ns,k) = etime(clock,ParticleStart);
                        else
                            local_fitness(ns) = ObjFuncHandle(LinPos(:,ns),ObjFuncData);
                        end
                    end
                end

            else

                if r.UseParallel
                    if r.ParticleElapsedTime
                        parfor ns = 1:r.NumSwp
                            ParticleStart = clock;
                            local_fitness(ns) = ObjFuncHandle(pos(:,ns),ObjFuncData); 
                            SwpElapsedTime(ns,k) = etime(clock,ParticleStart);
                        end
                    else
                        parfor ns = 1:r.NumSwp
                            local_fitness(ns) = ObjFuncHandle(pos(:,ns),ObjFuncData); 
                        end
                    end

                else
                    for ns = 1:r.NumSwp
                        if r.ParticleElapsedTime
                            ParticleStart = clock;
                            local_fitness(ns) = ObjFuncHandle(pos(:,ns),ObjFuncData); 
                            SwpElapsedTime(ns,k) = etime(clock,ParticleStart);
                        else
                            local_fitness(ns) = ObjFuncHandle(pos(:,ns),ObjFuncData);
                        end
                    end
                end
            end

        % Adjust the inertia coefficient
        if count <= r.wcount_min 
            w = min(w*r.wfactor_up,r.wmax);  % Must be less than wmax  
        end

        if count >= r.wcount_max
            w = max(w/r.wfactor_dn,r.wmin); % Must be greater than wmin
        end

        % Update the local bests and their fitness
        ind             = find(local_fitness < lf_best);
        lf_best(ind)    = local_fitness(ind);
        lp_best(:,ind)  = pos(:,ind);

        % Update global best and its fitness
        [min_l_fitness,index] = min(lf_best);

        if min_l_fitness < gf_best
            gf_best = min_l_fitness;  % Update global best fitness
            gp_best = pos(:,index);   % Update the global best param vector
            count   = max(count-1,0); % Decrease count by 1 (>= 0)
            l       = r.lmin;         % Decrease l to the lmin value

            % Save interim fit parameters to a mat file
            if ~isempty(r.FN_IntState)
                save(r.FN_IntState,'pos','vel','k','count','l','w','lf_best','gf_best','lp_best','gp_best','gf_bestHis','LogTransIdx','SwpElapsedTime');
            end

        else
            count   = max(count+1,0);% ensure count is always >= 0
            l       = min(l+r.lmin,r.NumSwp); % ensure l is always <= NumSwp
        end

        % Save current gf_best to a vector
        gf_bestHis(k) = gf_best;

        % Create a randomly selected neighborhood for each particle and
        % select the neighborhood particle with the best fitness.
        %
        %  Each local neighborhood is size "l" and contains "l"-1
        %  unique particles and the particle itself.
        %
        % TODO: collapse this loop if possible

        best_neigh_idx = zeros(r.NumSwp,1);
        for np = 1:r.NumSwp
            LocalNeigh = SwmNeigh_noSelf(randperm(r.NumSwp-1,l-1),np); 
            [LocalMinFit,LocalMinIdx] = min(lf_best(LocalNeigh));
            if LocalMinFit < lf_best(np)
                best_neigh_idx(np)  = LocalNeigh(LocalMinIdx);
            else 
                best_neigh_idx(np) = np;
            end
        end

        % Update particle velocities and positions
        vel = w*vel + r.c1*rand(NumPar,r.NumSwp).*(lp_best-pos)...
                    + r.c2*rand(NumPar,r.NumSwp).*(pos(:,best_neigh_idx)-pos);
        pos = pos + vel;  

            % Limit the new positions to values within the param bounds
            for np = 1:NumPar
                ind         = find(pos(np,:) <=r.LB(np));
                pos(np,ind) = r.LB(np);
                vel(np,ind) = 0;

                ind         = find(pos(np,:) >=r.UB(np));
                pos(np,ind) = r.UB(np);
                vel(np,ind) = 0;
            end

        % Evaluate early exit criteria

            % Break if the count is greater than the allowed number of
            % non-improving iterations (FitCount)
            if count > r.FitCount

                if r.Verbose
                     disp(['Exiting PSO loop - ' num2str(r.FitCount) ' count limit reached.']);
                end
                BreakFlag  = true;
                ExitReason = 'FitCount Reached';
                break

            end

        % Display status on screen
        if r.Verbose

            % Display
            if k == 1 || ~mod(k,25) 
                fprintf('  Iter \t  Fitness \t Count \t    l \t     w  \n')
            end

            fprintf('%6d \t %8.4e %6d \t %5d %11.3f \n',k,gf_best,count,l,w);
        end
    end

    % Linearize any log transformed parameter values
    if UseLogTransform
       gp_best(LogTransIdx) = 10.^(gp_best(LogTransIdx)); 
    end

    % Print fitness on screen and final fitness, parameters and velocities
    % on files
    if BreakFlag == false;
        ExitReason = 'MaxIter Reached';
        if r.Verbose;
            fprintf('Reached maximum number of iterations.\n Final fitness value is:');   fprintf('%12.4f\n',gf_best);
        end
    end
    
    % Calculate elapsed time
    EndClock    = clock;
    ElapsedTime = etime(EndClock,StartClock);

    if r.Verbose
        fprintf('%s %7.f %s\n','Elapsed time =',ElapsedTime,'seconds');
    end

    if ~isempty(r.FN_OptState)

        % If the user was saving interim results, simply rename the interim
        % results MAT file to the OptimalResults filename.  Else, create a
        % new MAT file.
        if ~isempty(r.FN_IntState) && exist(r.FN_IntState,'file') == 2
            movefile(r.FN_IntState,r.FN_OptState);
            % save(r.FN_OptState,'ExitReason','ElapsedTime','-append');
        else
            % save(r.FN_OptState,'pos','vel','k','count','l','w','lf_best','gf_best','lp_best','gp_best','gf_bestHis','LogTransIdx','SwpElapsedTime','ExitReason','ElapsedTime');
        end

    elseif ~isempty(r.FN_IntState) && exist(r.FN_IntState,'file') == 2
        % This would catch cases where the interim results are saved but
        % not the optimum results.
        % save(r.FN_IntState,'ExitReason','ElapsedTime','-append');
    end

end