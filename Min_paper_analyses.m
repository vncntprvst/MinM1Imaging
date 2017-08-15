%define parameters
preAlignWindow=2000; % 2 seconds before alignment time
postAlignWindow=5000; % 5 seconds after alignment time
samplingRate=10; %10Hz
interval=1000/samplingRate;
durationArray=1:(preAlignWindow+postAlignWindow)/interval;

subjectList={'F73';'F97';'F98';'F112'};
numSubject=size(subjectList,1);

%load data or
[data,averages,stats]=Min_paper_load_data(preAlignWindow,postAlignWindow,samplingRate,subjectList);

% allocate
params=struct('preAlignWindow',preAlignWindow,'postAlignWindow',postAlignWindow,...
    'samplingRate',samplingRate,'interval',interval,'durationArray',durationArray,...
    'numSubject',numSubject);
% figOptions=struct('alignSpecs',{},'figureStyle',{},'legends',{});
figOptions.alignSpecs={preAlignWindow/interval+1;...
    postAlignWindow/interval;10/interval;'white'};

%% extract calcium events for each epoch (already done, this is just to show)
% Min_paper_extract_calciumEvents(data);

%% plot reaching movement duration
% figure;
% for mvttplot=1:4
%     subplot(4,1,mvttplot)
%     plot(cellfun(@(mvttimes) diff(mvttimes), {data(mvttplot).behavData.movementTime}))
% end

%% run stats on neuron responses
%(takes a while if using bootstrap. Change boostrap parameter 'naccu' to 200 for faster operation)
stats=Min_paper_population_stats(data,params,figOptions);
% compare proportions classified cells
% proportions=...
%     [sum(vertcat(foo(1).taskRelated.indices))/numel(vertcat(foo(1).taskRelated.indices)),...
% sum(vertcat(foo(2).taskRelated.indices))/numel(vertcat(foo(2).taskRelated.indices)),...
% sum(vertcat(foo(3).taskRelated.indices))/numel(vertcat(foo(3).taskRelated.indices)),...
% sum(vertcat(foo(4).taskRelated.indices))/numel(vertcat(foo(4).taskRelated.indices));...
% ...
% ...%compare to bootstrap stats (all classifications pooled!):
% sum(([stats(1).taskRelated.indices]))/numel(([stats(1).taskRelated.indices])),...
% sum(([stats(2).taskRelated.indices]))/numel(([stats(2).taskRelated.indices])),...
% sum(([stats(3).taskRelated.indices]))/numel(([stats(3).taskRelated.indices])),...
% sum(([stats(4).taskRelated.indices]))/numel(([stats(4).taskRelated.indices]))]
% 
% figure; bar(proportions'); legend({'signrank','bootstrap'})

%% categories (pie charts)
catPct=Min_paper_categorize_neurons(stats);

%% display population responses
averages=Min_paper_population_response_profiles(data,params,foo,figOptions); %fill in [] for stats if none available

%% is variability decreasing?
Min_paper_variability(data,stats,params,figOptions);

%% Figures for Saturday
% a. population timing across sessions
% done with Min_paper_population_response_profiles
% b. individual neurons are not behavior outcome specific. For individual neurons,
%     timing and behavior outcome are not correlated
% run the  %% test movement response vs outcome in Min_paper_population_stats
%     (stats saved in stats_on_caEvents_with_direction_classbysigTime.mat)
% -> find non-random association between response at reach time and outcome
outcomeRelated=[stats.outcomeRelated];
responseCorOutcome{1}=[outcomeRelated(1:3:end).indices];
responseCorOutcome{2}=[outcomeRelated(2:3:end).indices];
responseCorOutcome{3}=[outcomeRelated(3:3:end).indices];
figure;bar([sum(responseCorOutcome{1})/numel(responseCorOutcome{1}),...
    sum(responseCorOutcome{2})/numel(responseCorOutcome{2}),...
    sum(responseCorOutcome{3})/numel(responseCorOutcome{3})]);
set(gca,'xticklabel', {'Naive Session','Expert 1','Expert 2'},'ylim',[0 0.1],...
    'FontWeight','bold','FontSize',12);
ylabel('Percentage of cells associated to outcome')
colormap('bone');

% c. trial by trial variation both temporarily and in percentage of activated neurons,
%     if possible, could you please compare ‘s’ trials to all other trials?
% -> use Min_paper_variability








