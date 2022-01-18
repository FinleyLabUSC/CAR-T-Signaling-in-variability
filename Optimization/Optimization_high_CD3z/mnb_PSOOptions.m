%% mnb_PSOOptions 
% Get options struct for particle swarm optimization (PSO)


function options = mnb_PSOOptions()

    options = struct;

    % Objective Function options
    options.ObjFuncDir  = '';   % char, path to the directory containing the objfunc.
                                %   If empty, the current path (pwd) will be used.
    options.ObjFuncData = [];   % user-defined variable that is passed to the objfunc.
                                %   It can be of any type (struct, cell, etc) and is
                                %   meant to be used to pass the experimental data
                                %   and other flags to the objective function.

    % Initial Values for parameters
    %  mnb_PSOFit requires lower and upper bounds for parameter values to
    %  be supplied as function input arguments. The fit routine is then
    %  initialized with a random population between these bounds. The user
    %  can also supply an initial guess (vector of parameter values) that
    %  will be included as a member of the initial population.
    options.InitialGuess  = []; % vector of length #_of_parameters_being_fit

    % Log Transform parameter values
    %  mnb_PSOFit can internally log-transform certain or all parameter
    %  values to better search a broad parameter space. This is transparent
    %  to the user, as they supply linear bound values (LB & UB) and linear
    %  values are passed to the objective function.
    %
    % To enable log-transformation, define the following Param-Value
    % argument as a logical vector of length #Parameters. Indices set to
    % TRUE will be log-transformed while those set to FALSE will not.
    options.LogTransform = []; % if empty, no parameters will be log-transformed

    % PSO loop exit criteria

        % Exit Criterion 1: Stop fitting at a pre-determined number of
        % iterations. This is only invoked if the other early-exit criteria
        % (below) are not achieved.
        options.MaxIter     = 2500;  % # iterations

        % Exit Criterion 2: Stop fitting if the number of non-improving
        % iterations (captured in the "count" variable) reaches this value.
        options.FitCount    = 100; % # iterations

    % Swarm parameters

        % Number of swarm particles
        options.NumSwp = 20; % positive integer >=1

        % Randomly assigned neighborhoods
        options.lmin = 5;

        % Min and max inertia factors
        options.wmin = 0.1;   
        options.wmax = 1.1;   

        % Min and max values of counter for increasing and decreasing inertia
        options.wcount_min = 4;   % integer >= 1
        options.wcount_max = 7;   % integer > wcount_max

        % Scaling values for the adaptive inertia factor. The intertia
        % factor will be increased (w*wfactor_up) when count > wcount_min
        % and decreased (w/wfactor_dn) when count > wcount_max
        options.wfactor_up = 1.5; % numeric scalar > 0
        options.wfactor_dn = 1.5; % numeric scalar > 0

        % Self-recognition and Social component coefficients
        options.c1 = 1.49;   
        options.c2 = 1.49;   

        % Status display and file exports 

        % Display status updates for each iteration
        options.Verbose = false; 

        % File Name for final output file (MAT) : fitness and best params
        options.FN_OptState = 'PSO Fit - Final results.mat';

        % File Name for full interim status (MAT) : all PSO fit loop
        % variables necessary to restart a PSO routine from the same state.
        options.FN_IntState = '';          % 'PSO Fit - interim fit state.mat'

        % To restart from a prior fit, supply the file name for the MAT
        % file containing the saved state.
        options.FN_IntState_restart  = ''; % 'Interim fit state from last fit.mat'

        % Extend a simulation that hit the MaxMaxIter number
        options.ExtendMaxIter = false;

    % Other options

        % Randstream initialization
        options.InitRandStream = true; % initialize the randstream

        % Use the parallel computing toolbox : call the objective function
        % within a parfor loop (true) or a regular for loop (false)
        options.UseParallel = true;

        % Save time elapsed by each particle
        options.ParticleElapsedTime = false;
