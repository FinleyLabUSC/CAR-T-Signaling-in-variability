%% mnb_PSOFitSummary
% Generate diagnostic plots of multiple, independent particle swarm
% optimization (PSO) fits using the .MAT files saved from each fit. This
% function will load all .MAT files contained within the directory
% specified, so be careful to remove any extraneous .MAT files beforehand.
%
% Usage: [All_FitParams,All_FitVals,All_FitValsHist] = mnb_PSOFitSummary(FitDir,varargin)
%

%% Required Input Arguments (in order)
%
% * FitDir : char : path to the directory containing the .MAT files. 
%
%% Optional Input Arguments (passed via Param-Value pairs)
%
% * LB : numeric vector : lower bounds for estimated parameter values. 
%
% * UB : numeric vector : upper bounds for estimated parameter values. 
%
% * FitParamNames : cell string : names of estimated parameters.
%
% * FitValThreshold : numeric scalar : when plotting the coefficient of
%   variation for each parameter, also include a plot of the subset of fits
%   with a fitness value less than this threshold.
%
%% Function Output Arguments
%
% * All_FitParams : numeric array : concatenation of parameter values
%       loaded from all .MAT files. The size of the array is the number of
%       parameters x the number of .MAT files.
%
% * All_FitVals : numeric vector : concatenation of fitness values loaded
%       from all .MAT files. The length is equal to the number of .MAT
%       files.
%
% * All_FitValsHist : numeric array : concatenation of fitness value
%       trajectories (fitness value at each iteration of the PSO fit
%       routine) loaded from all .MAT files. The size of the array is the
%       max number of iterations x the number of .MAT files.
%
%% Examples
%
%  % Generate plots for fitness values, time to convergence, and param CVs.
%  FitDir = '\some\path\to\a\folder\with\fit\results';
%  mnb_PSOFitSummary(FitDir);
%
%  % Do not generate plots. Only return the output arguments.
%  FitDir = '\some\path\to\a\folder\with\fit\results';
%  mnb_PSOFitSummary(FitDir,'GenPlots',false);
%
%% Dependencies
%
% * Statistics toolbox for nanstd and nanmean functions
%
%% Limitations
%
% * Need to update the parameter plots to include upper and lower bounds.
% This will help identify instances where parameter bounds should be
% adjusted.
%
% * It would be useful to identify and plot parameter ranges for the global
% versus local minima. Or at the least, allow the user to define fitness
% value ranges for these minima and have subplots generated for each.
%
%% Version Information
% * Last changed date: $Date:$
% * Last changed by: $Author:$
% * Last changed revision: $Revision:$
%
% Copyright 2013-2015 Merrimack Pharmaceuticals, Inc.

function [All_FitParams,All_FitVals,All_FitValsHist] = mnb_PSOFitSummary(FitDir,varargin)

%-- Validate and parse the input
    % Configure the Parser
    p = inputParser();
    p.FunctionName  = 'MNB_PSOFITSUMMARY';
    p.CaseSensitive = false;
    p.StructExpand  = true;
    p.KeepUnmatched = true; %ignore unknown input arguments

    % Validate Required Input Arguments (must be passed in this order)
    p.addRequired('FitDir',@(x) ischar(x) && exist(x, 'dir'));

    % Validate Optional Input Arguments (use Param-Value pairs)
    p.addParameter('LB',[],@(x) isnumeric(x) && ~isempty(x) && (isrow(x) || iscolumn(x)));
    p.addParameter('UB',[],@(x) isnumeric(x) && ~isempty(x) && (isrow(x) || iscolumn(x)));
    p.addParameter('FitParamNames',[],@(x) ~isempty(x) && isvector(x) && iscellstr(x));
    p.addParameter('FitValThreshold',[],@(x) isnumeric(x) && isscalar(x));
    p.addParameter('GenPlots',true,@islogical);

    % Parse the input arguments and push to the struct 'r'
    p.parse(FitDir,varargin{:});  
    r = p.Results();

    % Post-processing of parsed input arguments (mostly error checking)

    if ~isempty(r.LB) && ~isempty(r.UB)
        
        GenParamPlots = true;
        
        % Convert LB, UB, and FitParamNames to column (if row)
        if isrow(r.LB); r.LB = r.LB'; end
        if isrow(r.UB); r.UB = r.UB'; end
        
        if ~isempty(r.FitParamNames) && isrow(r.FitParamNames)
            r.FitParamNames = r.FitParamNames';
        end
        
    else
        GenParamPlots = false;
    end

%-- Build a list of .MAT files in the FitDir directory
    
    MatFiles = dir(fullfile(r.FitDir,'*.mat'));
    NumFiles = length(MatFiles);
    
%-- Load the first .MAT file in order to determine the number of Max
%       Iterations and the number of fit parameters.
    
    try
        load(fullfile(r.FitDir,MatFiles(1).name),'gf_best','gf_bestHis','ElapsedTime','gp_best');

    catch % something wrong with mat file or variables missing
        error([p.FunctionName ':1stFileLoadError'],['Failed to load file "' MatFiles(1).name '"']);
    end

    MaxIter     = length(gf_bestHis);
    NumFitParams= length(gp_best);

        if ~isempty(r.LB) && (NumFitParams ~= length(r.LB))
            error([p.FunctionName ':NumParams'],'The # of parameters specified in "LB" does not match "gp_best" in the .MAT file.');
        end

    clear gf_best gf_bestHis ElapsedTime gp_best;
    
%-- Load each file and compile results   
    All_FitVals     = nan(NumFiles,1);
    All_FitValsHist = nan(MaxIter,NumFiles);
    All_TimeElapsed = nan(NumFiles,1);
    All_FitParams   = nan(NumFitParams,NumFiles);

    for nf = 1:NumFiles

        % Load File
        try
            load(fullfile(r.FitDir,MatFiles(nf).name),'gf_best','gf_bestHis','ElapsedTime','gp_best','LogTransIdx');

        catch % something wrong with mat file or variables missing
            disp(['Fit File #' num2str(nf) ' load error']);
            break;
        end

        % save file-specific variables to arrays
        All_FitVals(nf)         = gf_best;
        All_FitValsHist(:,nf)   = gf_bestHis;
        All_TimeElapsed(nf)     = ElapsedTime;
        
            % Linearize any log-transformed parameter values
            gp_best(LogTransIdx) = 10.^gp_best(LogTransIdx);
            All_FitParams(:,nf)  = gp_best;
        
        % clear variables (setup for next load)
        clear gf_best gf_bestHis ElapsedTime gp_best LogTransIdx

    end
  
%-- Plot best fitness values and elapsed time for each separate fit
    if r.GenPlots
        figure;
        [gf_bestSorted,gf_bestSortOrder] = sort(All_FitVals);

        subplot(3,1,1);
        plot(All_FitVals(gf_bestSortOrder),'o');
        title('Sorted Final Fitness Values (All Runs)');
        ylabel('Fitness Value');
        xlabel('Fit Run (sorted)');

        % Plot Time Elapsed for each fit run
        subplot(3,1,2);
        plot(All_TimeElapsed(gf_bestSortOrder)/60,'o');
        title('Time to Convergence sorted by FitVal (All Runs))');
        ylabel('Time Elapsed (min)');
        xlabel('Fit Run (sorted)');

        % Plot gf_bestHis for each fit run
        subplot(3,1,3)
        for nf = 1:NumFiles
            hold on; plot((1:MaxIter),All_FitValsHist(:,nf),'Color',[0.6 0.6 0.6]);
        end

        % Plot the median of all values
        hold on; plot((1:MaxIter),nanmedian(All_FitValsHist,2),'r');

        title('Fitness Values by iteration # (All Runs)');
        ylabel('Fitness Value');
        xlabel('Iteration #');

            % Set the limits on the X-axis
            NumFits = size(All_FitValsHist,2);

            LastNaNIdx = nan(NumFits,1);

            for nf = 1:NumFits
                NaNIdx = find(isnan(All_FitValsHist(:,nf)),1,'first');

                if isempty(NaNIdx)
                    LastNaNIdx(nf) = MaxIter;
                else
                    LastNaNIdx(nf) = NaNIdx;
                end
            end

            XMax = max(LastNaNIdx);

        set(gca,'XLim',[0 XMax]);


    %-- Evaluate variability in parameter values
        fh2 =  figure; 

        % How many plots?
        if ~isempty(r.FitValThreshold)
            NumSubplots = 2;
        else
            NumSubplots = 1;
        end

        % Parameter CV -- all fits
        ParamCV = nanstd(All_FitParams,1,2) ./nanmean(All_FitParams,2);
        [SortedParamCV,CVorder] = sort(ParamCV);

        subplot(1,NumSubplots,1);
        barh(SortedParamCV);
        title('CV value for each Parameter');
        set(gca,'YTick',(1:1:length(ParamCV)));
        xlabel('CV Value -- nanstd/nanmean');

        if ~isempty(r.FitParamNames) && length(r.FitParamNames) == length(ParamCV)
            set(gca,'YTickLabel',r.FitParamNames(CVorder));
        end

        % Parameter CV -- subset of fits with fitness < threshold
        if ~isempty(r.FitValThreshold)
            subplot(1,NumSubplots,2);

            FitValSubset = All_FitVals <= r.FitValThreshold;

            if any(FitValSubset)
                ParamCV = nanstd(All_FitParams(:,FitValSubset),1,2) ./nanmean(All_FitParams(:,FitValSubset),2);
                [SortedParamCV,CVorder] = sort(ParamCV);

                subplot(1,NumSubplots,2);
                barh(SortedParamCV);
                title(['CV value for each Parameter, Fitness <' num2str(r.FitValThreshold)]);
                set(gca,'YTick',(1:1:length(ParamCV)));
                xlabel('CV Value -- nanstd/nanmean');

                if ~isempty(r.FitParamNames) && length(r.FitParamNames) == length(ParamCV)
                    set(gca,'YTickLabel',r.FitParamNames(CVorder));
                end
            end
        end

    end
end