% Written by Keo 20200524
% including 'heading angles'

clc;clear;
% =========================================== ¡ï¡ï¡ï¡ï¡ï =========================================
% === PSTH, two monkeys ===
FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20191010_Messi_and_Fara_FEF_dt_PSTH_bin10_CorrectTrial\';% binSize = 10; correct trials, Messi & Fara

cond_num = 3;% 3 - no dt;5 - include dt
             % 1:ves, 2:vis, 3:comb0, 4:comb-500, 5:comb-250

LEFT = 1;
RIGHT = 2;

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath,'*.mat'));
File_work_count = 0;
list_1 = load([FolderPath Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_1,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name{File_work_count} = Files(ct).name;
        list{File_work_count} = load([FolderPath Files(ct).name]);
    end
end
File_normal_total = File_work_count;
display(['Total = ',num2str(File_normal_total)])

% files to analyze
file_to_analyze = 1:File_normal_total;
% file_to_analyze = [98];
display(['To analyze = ',num2str(length(file_to_analyze))])

if length(file_to_analyze) >= 2
    isGroup = 1;
elseif length(file_to_analyze) == 1
    isGroup = 0;
end

% 0.2 colors, markers
symbol_color{1} = [0 0 1]; % blue - vestibular only
symbol_color{2} = [1 0 0]; % red - visual only
symbol_color{3} = [0 1 0];% green - combined in 'dt=0'
symbol_color{4} = [0 1 1];%[51 204 204]/255;% cyan
symbol_color{5} = [1 0 1];%[255 51 204]/255;% magenta
symbol_color{6} = [1 1 0];%[255 134 40]/255;% orange

symbol{1} = '-';
symbol{2} = '--';

line_width_axis = 1;% axis
line_width_curve = 1;% curve
line_width_vel = 0.8;% velocity profile

font_size = 15;

legend_str_dt = {' Ves ' '  Vis ' '     0 ' '-500' '-250'};
legend_str_PN = {'Pref ¡ª' 'Null ---'};

%% find the max trial_num
max_trial_num = 0;
j = 1; % 1 = Stim On
for f = 1:length(file_to_analyze)
    for cond = 1:cond_num % 1:ves, 2:vis, 3:comb0, 4:comb-500, 5:comb-250
        for PN = 1:2 % 1:pref, 2:null
            cond_idx = list{1,f}.result.stim_type_per_trial' == cond;
            if list{1,f}.result.PREF == 1
                PN_idx = list{1,f}.result.choice_per_trial == PN;
            else
                PN_idx = ~(list{1,f}.result.choice_per_trial == PN);
            end
            trial_num = sum(cond_idx & PN_idx  == 1);
            if trial_num > max_trial_num
                max_trial_num = trial_num;
            end
        end
    end
end

%%
ts = list{1,f}.result.rate_ts{1, 1};
ts_num = size(ts,2);
j = 1; % 1 = Stim On
tic
clear firingRates
for f = 1:length(file_to_analyze)
    for cond = 1:cond_num % 1:ves, 2:vis, 3:comb0, 4:comb-500, 5:comb-250
        for PN = 1:2 % 1-pref, 2-null
            cond_idx = list{1,f}.result.stim_type_per_trial' == cond;
            if list{1,f}.result.PREF == 1
                PN_idx = list{1,f}.result.choice_per_trial == PN;
            else
                PN_idx = ~(list{1,f}.result.choice_per_trial == PN);
            end
            trial_num = sum(cond_idx & PN_idx == 1);
            
            firingRates_padding = [list{1,f}.result.spike_hist{1,j}(cond_idx & PN_idx,:); nan(max_trial_num - trial_num, ts_num)]';
            firingRates(f, cond, PN, 1:ts_num, 1:max_trial_num) = firingRates_padding;
            
        end
    end
end
toc
firingRatesAverage = nanmean(firingRates,5);

%% plot average PSTH
COLORS = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 0.7 0.7; 0.7 1 0.7];
set(figure(123),'color','white');clf;hold on
for st = 1:size(firingRatesAverage,2)
    for PN = 1:size(firingRatesAverage,3)
        PSTH = mean(squeeze(firingRatesAverage(:, st, PN, :)),1);
        if PN == 1
            plot(PSTH, 'color', COLORS(st,:), 'LineWidth', 2,'linestyle','-');
        else
            plot(PSTH, 'color', COLORS(st,:), 'LineWidth', 2,'linestyle','--');
        end
    end
end
%%
% % N is the number of neurons
% % S is the number of stimuli type conditions
% % D is the number of decisions (D=2)
% % T is the number of time-points (note that all the trials should have the same length in time!)
% 
% N = length(Files);   % number of neurons
% S = 5;               % number of stimuli type
% D = 2;               % number of decisions
% T = length(ts);      % number of time points
% E = max_trial_num;   % maximal number of trial repetitions

%% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus type', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
% margColours = [1 0 0; 0 1 0; 0 0 1; 1 1 0];

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
[~,timeEvents] = min(abs(ts));
timeEvents = [0 2000];

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', ts,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
    'combinedParams', combinedParams);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', ts,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);


%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations
% in a .mat file with a given name. Once computed, you can simply load
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', ts,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', ts,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: decoding

decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};

accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], [], [], decodingClasses)

accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', ts,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif);
