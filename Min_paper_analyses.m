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

%% run stats on neuron responses
%(takes a while. Change boostrap parameter 'naccu' to 200 for faster operation)
stats=Min_paper_population_stats(data,params,figOptions);

%% categories (pie charts)
catPct=Min_paper_categorize_neurons(stats);

%% display population responses
averages=Min_paper_population_response_profiles(data,params,stats,figOptions); %fill in [] for stats if none available

%% is variability decreasing? 
Min_paper_variability(data,stats,params);

%% Figures for Saturday
% a. population timing across sessions

% b. individual neurons are not behavior outcome specific. For individual neurons,
%     timing and behavior outcome are not correlated
% c. trial by trail variation both temporarily and in percentage of activated neurons,
%     if possible, could you please compare ‘s’ trails to all other trials?

% does learning result in sparser representation? 



  
    
   


