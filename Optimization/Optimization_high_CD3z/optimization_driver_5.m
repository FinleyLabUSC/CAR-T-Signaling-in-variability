 userpath='/project/sfinley_292/vtseruny/NFkB/Optimization_high_CD28';
 MPCTprofile=parallel.cluster.Local;
 MPCTprofile.NumWorkers=21;
 MPCTprofile.JobStorageLocation='/project/sfinley_292/vtseruny/Optimization_high_CD28';
 parpool(MPCTprofile, 20);

tic

% Load parameter bounds and experimental data
load('Bounds.mat');
LB = Bounds(:,1);
UB = Bounds(:,2); 
gp_aggregator=nan(25,length(LB));
gf_aggregator=nan(25,1);
load('ObjFuncData.mat');


% Configure the PSO options
PSOopts = mnb_PSOOptions();

% PSO Configuration 
PSOopts.ObjFuncData = ObjFuncData;
PSOopts.LogTransform = true(length(LB),1);
PSOopts.NumSwp  = 20;  % Number of swarm particles = number of the CPUs
PSOopts.lmin    = 3;
PSOopts.MaxIter     = 1000;
PSOopts.FitCount    = 50;   % Exit if PSO count reaches this value
PSOopts.Verbose         = true; % disable status display
PSOopts.FN_IntState     = ''; % 'EarlierIntState.mat';
PSOopts.UseParallel  =  true;

% Run PSO for a set number of times
for i=1:11

    ObjFuncName = 'objective_function_5';
    [gp_best,gf_best] = mnb_PSOFit(ObjFuncName,LB,UB,PSOopts);
    display(gp_best)
    gp_aggregator(i,:)=gp_best;
    gf_aggregator(i)=gf_best;

    save('PSOresults_5.mat','gp_aggregator','gf_aggregator')

end

toc
