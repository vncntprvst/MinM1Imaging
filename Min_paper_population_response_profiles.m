function averages=Min_paper_population_response_profiles(data,params,stats,figOptions)

durationArray=params.durationArray;
preAlignWindow=params.preAlignWindow;
postAlignWindow=params.postAlignWindow;
interval=params.interval;
numSubject=params.numSubject;

%% sanity check: draw first spike from first cell, 
% % from full raw trace vs from the trial epoch's raw trace
% figure;
% spikeOneTrial=find(cellfun(@(x) ~isempty(x{1}), {data(subjectNum).ciData.spikes}),1);
% % define time window
%     targetTime(1)=data(subjectNum).ciData(spikeOneTrial).spikes{1, 1}(1)-10999;
%     targetTime(2)=data(subjectNum).ciData(spikeOneTrial).spikes{1, 1}(2)+20001;
% % find corresponding data indices
%     windowIdx(1)=find(data(subjectNum).rawTraces{1,1}.Time==targetTime(1)); % assuming first spike is in first session.  
%     windowIdx(2)=find(data(subjectNum).rawTraces{1,1}.Time==targetTime(2)); % that might not always be true
% %plot zoomed out view
% subplot(2,1,1)
% plot(data(subjectNum).rawTraces{1,1}.Data(windowIdx(1):windowIdx(2),1));
% hold on 
% % overlap spike trace
% spikeTimeArray=windowIdx(1)+110:windowIdx(2)-200;
% plot(111:windowIdx(2)-windowIdx(1)-199,data(subjectNum).rawTraces{1,1}.Data(spikeTimeArray,1));
% yLims=get(gca,'ylim');
% % now get trial traces
% targetTime=data(subjectNum).ciData(spikeOneTrial).trialWindowIndex;
% windowIdx(1)=find(data(subjectNum).rawTraces{1,1}.Time>=targetTime(1),1);
% windowIdx(2)=find(data(subjectNum).rawTraces{1,1}.Time<=targetTime(2),1,'last');
% % trial plot
% subplot(2,1,2) 
% % plot([zeros(100,1);data(subjectNum).ciData(spikeOneTrial).rawTraceEpochs.Data(:,1);zeros(200,1)])%zero-padded trace for rough adjustment
% plot(data(subjectNum).rawTraces{1,1}.Data(windowIdx(1):windowIdx(2),1));
% hold on
% % overlap spike trace
% spikeTimeArray=find(data(subjectNum).ciData(spikeOneTrial).rawTraceEpochs.Time>=...
% data(subjectNum).ciData(spikeOneTrial).spikes{1, 1}(1),1):...
% find(data(subjectNum).ciData(spikeOneTrial).rawTraceEpochs.Time<=...
% data(subjectNum).ciData(spikeOneTrial).spikes{1, 1}(2),1,'last');
% plot(spikeTimeArray,data(subjectNum).ciData(spikeOneTrial).rawTraceEpochs.Data(spikeTimeArray,1)) 
% % set(gca,'ylim',yLims);
 
%% display activity rasters for one particular cell
cellNum=36;
subjectNum=4;
% cell index (for cells that could be tracked across sessions)
cellIndex=table2array(data(subjectNum).cellIDs(1).cellIndex);

% find bad / no spike / too short trials
badTrials=cellfun(@(spikeCell) size(spikeCell,2)==1,...
    {data(subjectNum).ciData.spikes})'; %likely interrupted trial at the end of a session
noSpikeTrials=~cellfun(@(spikeCell) sum(cellfun('isempty',spikeCell))<size(spikeCell,2),...
    {data(subjectNum).ciData.spikes});
sampleLength=cell2mat(cellfun(@(trialResponse) get(trialResponse,'Length'),...
    {data(subjectNum).ciData.rawTraceEpochs}','UniformOutput', false));
shortTrials=sampleLength<max(durationArray); %to exclude trials that are too short 

% allocate
ciTraces=NaN(size(sampleLength,1),numel(durationArray));
ciSpikes=cell(size(sampleLength,1),1);

% session index
% e.g., keep expert sessions only
sessionDays=unique([data(subjectNum).behavData.session]);
% allSessionIdx=[data(subjectNum).behavData.session]'== sessionDays(2) |...
%      [data(subjectNum).behavData.session]'== sessionDays(3);
keepSessions=[2,3];

for sessionNum=keepSessions(1):keepSessions(end)
    sessionIdx=[data(subjectNum).behavData.session]'==sessionDays(sessionNum);
    trialIdx=~(badTrials | shortTrials | ~sessionIdx);
    cellId=cellIndex(cellIndex(:,2)==cellNum,sessionNum); 
    
    sTraces=data(subjectNum).rawTraces{sessionNum, 1}.Data;
    trialTimes={data(subjectNum).behavData(trialIdx).movementTime};
    
    ciTraces(trialIdx,:)=cell2mat(cellfun(@(mvtimes) zscore(sTraces(...
        durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellId)') ,...
        trialTimes,'UniformOutput', false)');
    
%     foo=cell2mat(cellfun(@(trialResponse) zscore(trialResponse.Data(durationArray,cellId)') ,...
%         {data(subjectNum).ciData(trialIdx).rawTraceEpochs}','UniformOutput', false));
%     figure; hold on; plot(foo(end,:)); plot(ciTraces(find(trialIdx,1,'last'),:))
    
    % spikes, referenced to defined time window
%     ciSpikes(trialIdx,:)=cellfun(@(trialSpikes,trialTimeWindows)...
%         (trialSpikes{cellId}-trialTimeWindows(1))/interval,...
%         {data(subjectNum).ciData(trialIdx).spikes},{data(subjectNum).ciData(trialIdx).trialWindowIndex},...
%         'UniformOutput', false);
    % spikes, referenced to actual epoch's start time (not much of a difference mind you)
    ciSpikes(trialIdx,:)=cellfun(@(trialSpikes,mvtimes)...
        (trialSpikes{cellId}-(mvtimes(1)-preAlignWindow))/interval,...
        {data(subjectNum).ciData(trialIdx).spikes},trialTimes,...
        'UniformOutput', false);
end

% plot one single trial 
% trialResponse={data(subjectNum).ciData(trialIdx).rawTraceEpochs}';
% trialSpikes={data(subjectNum).ciData(trialIdx).spikes};
% trialTimeWindows={data(subjectNum).ciData(trialIdx).trialWindowIndex};
% trial=1;
% trace=zscore(trialResponse{trial}.Data(durationArray,cellId)');
% figure;plot(trace);
% spikeTimes=(trialSpikes{trial}{cellId}-trialResponse{trial}.Time(1))/interval;
% patch(sort(repmat(spikeTimes,1,2)), [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
%     [0 0 0 0],[0.4 0.3 0.7],'EdgeColor','none','FaceAlpha',0.3);
% spikeTimes=(trialSpikes{trial}{cellId}-trialTimeWindows{trial}(1))/interval;
% patch(sort(repmat(spikeTimes,1,2)), [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
%     [0 0 0 0],[0.7 0.3 0.4],'EdgeColor','none','FaceAlpha',0.3);

% remove spurious rows (bad trials, not in session, etc)
ciSpikes=ciSpikes(~isnan(sum(ciTraces,2)),:);
ciTraces=ciTraces(~isnan(sum(ciTraces,2)),:);

% plot rasters
figOptions.figureStyle={'parula';1:10:size(ciTraces,2)+1;-preAlignWindow/1000:postAlignWindow/1000};
figOptions.legends={'Time (s)'; 'Trials'; {['Neuron # ' num2str(cellNum)],...
        'calcium response aligned to movement initiation'}};
% plot_rasters_psth({ciTraces,ciSpikes},figOptions);

%%%%%%%%%%%%%%%%%%%%%%

%% Population plots %%
%%%%%%%%%%%%%%%%%%%%%%

%% Trial averages of all classified neurons
% and find task-related neurons 

% Li et al 2015: use Kruskal–Wallis across multiple 0.5s bins (five image samples per bin) 
% null hypothesis: all time bins have equal fluorescence. Reject at P < 0.01 

% ?? : Wilcoxon signed-rank test on average of ±2 frames 

% Peters 2014: The amount of activity during movement
% was calculated for each neuron as the mean value of the activity event trace during
% movement epochs (defined as described above, and extending individual epochs
% by 5 image frames before and after each movement). The movement trace was then
% shuffled (10,000 times) such that complete movement epochs were kept intact but
% their position in the trace and relation to each other was randomized.Ameasure of
% activity during these shuffled movement epochs was calculated in each shuffle as
% above. The neuron was classified as movement-related if the real value was higher
% than the 0.5 percentile value of the shuffled values.

%Define task epoch 5 frames before movement and 3 after (8 total), grouped in 4 bins, then same as Li et al. 2015

% [averages.trialAverages,data.taskRelated]=deal(cell(3,2));%3 sessions: 1 naive, 2 experts 
averages=struct(data); fldName=fieldnames(averages);
averages=rmfield(averages,fldName(~cellfun(@(x) contains(x,'subject'), fldName)));
[averages.trialAverages]=deal(cell(3,2));%3 sessions: 1 naive, 2 experts 

for subjectNum=1:numSubject
    % find bad / no spike / too short trials
    badTrials=cellfun(@(spikeCell) size(spikeCell,2)==1,...
        {data(subjectNum).ciData.spikes})'; %likely interrupted trial at the end of a session
    sampleLength=cell2mat(cellfun(@(trialResponse) get(trialResponse,'Length'),...
        {data(subjectNum).ciData.rawTraceEpochs}','UniformOutput', false));
    shortTrials=sampleLength<max(durationArray); %to exclude trials that are too short
    sessionDays=unique([data(subjectNum).behavData.session]);
    for sessionNum=1:numel(sessionDays)
        sessionIdx=[data(subjectNum).behavData.session]'==sessionDays(sessionNum);
        trialIdx=~(badTrials | shortTrials | ~sessionIdx);
        
        %keep indices for each outcomes %m:missed f:failed n:no seed d:droped s:success
        trialOutcomes={data(subjectNum).behavData(trialIdx).outcome}; 
%         [outcomeIdx.missed,outcomeIdx.failed,outcomeIdx.noseed,outcomeIdx.droped,outcomeIdx.success]=deal(trialIdx);
        outcomeIdx.missed=cellfun(@(trialType) contains(trialType,'m'), trialOutcomes);
        outcomeIdx.failed=cellfun(@(trialType) contains(trialType,'f'), trialOutcomes);
        outcomeIdx.noseed=cellfun(@(trialType) contains(trialType,'n'), trialOutcomes);
        outcomeIdx.droped=cellfun(@(trialType) contains(trialType,'d'), trialOutcomes);
        outcomeIdx.success=cellfun(@(trialType) contains(trialType,'s'), trialOutcomes);

        %
        movementTimes={data(subjectNum).behavData(trialIdx).movementTime};
        numCells=size(data(subjectNum).ciData(find(trialIdx,1)).spikes,2);
        
        % extract calcium traces and events for that session's trials
        caRawTrace=data(subjectNum).rawTraces{sessionNum,1}.Data;
        caEventTrace=data(subjectNum).calciumEvents{sessionNum};
%         trialTimeArrays=cellfun(@(trialArray)(trialArray.Time(1:max(durationArray))-1)/100,...
%         {data(subjectNum).ciData(trialIdx).rawTraceEpochs},'UniformOutput', false);  
    
        % run analysis with zscore on whole trace, then compare when using
        % baseline zscore normalization (using trace before movement initation)
        
        [ciTraces,ciTraces_noz,caEvents,caEvents_noz]=deal(NaN(numCells,numel(durationArray)));
%         data(subjectNum).taskRelated{sessionNum,1}=NaN(numCells,3); %col1: is task related 
                            %/ col2: first time of significant difference / col3: duration 
%         data(subjectNum).taskRelated{sessionNum,2}=cell(numCells,1); % col4: classification

        for cellNum=1:numCells
            % average raw traces
            ciTraces(cellNum,:)=zscore(mean(cell2mat(cellfun(@(mvtimes) caRawTrace(...
                durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
                movementTimes,'UniformOutput', false)')));
            ciTraces_noz(cellNum,:)=mean(cell2mat(cellfun(@(mvtimes) caRawTrace(...
                durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
                movementTimes,'UniformOutput', false)'));
            %             ciTraces(cellNum,:)=zscore(mean(cell2mat(cellfun(@(trialResponse) trialResponse.Data(durationArray,cellNum)' ,...
            %                 {data(subjectNum).ciData(trialIdx).rawTraceEpochs}','UniformOutput', false))));
            
            %averages events
            caEvents(cellNum,:)=zscore(mean(cell2mat(cellfun(@(mvtimes) ...
                caEventTrace(durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum),...
                movementTimes,'UniformOutput', false)),2));
            caEvents_noz(cellNum,:)=mean(cell2mat(cellfun(@(mvtimes) ...
                caEventTrace(durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum),...
                movementTimes,'UniformOutput', false)),2);
            %figure; hold on; plot(caEvents(cellNum,:)); plot(ciTraces(cellNum,:))
        end
%         figure; hold on;
%         plot(ciTraces(62,:)); plot(cEvents(62,:));
        averages(subjectNum).trialAverages{sessionNum,1}=ciTraces;
        averages(subjectNum).trialAverages{sessionNum,2}=caEvents;
        averages(subjectNum).trialAverages{sessionNum,3}=ciTraces_noz;
        averages(subjectNum).trialAverages{sessionNum,4}=caEvents_noz;
    end
end

%% clear variables
clearvars -except data averages figOptions durationArray preAlignWindow postAlignWindow stats

%% plot each session, with traces ordered by the timing of their peak activity
% sessionTraces{1}=[averages.trialAverages];sessionTraces{1}=sessionTraces{1}(1,1:4:end)';
% sessionTraces{1}=vertcat(sessionTraces{1}{:});
% % sessionTaskRelatedCells{1}=[data.taskRelated];sessionTaskRelatedCells{1}=sessionTaskRelatedCells{1}(1,1:2:end)';
% % sessionTaskRelatedCells{1}=vertcat(sessionTaskRelatedCells{1}{:});
% [~,peakIdx]=max(sessionTraces{1},[],2);
% [~,peakOrder]=sort(peakIdx);
% % sessionTaskRelatedCells{1}=sessionTaskRelatedCells{1}(peakOrder,1);
% sessionTraces{1}=sessionTraces{1}(peakOrder,:);
% 
% figOptions.figureStyle={'parula';1:10:size(durationArray,2)+1;-preAlignWindow/1000:postAlignWindow/1000};
% figOptions.legends={'Time (s)'; 'Cells'; {'Session 1 (naive) ',...
%         'calcium response aligned to movement initiation'}};
% args{1,1}={sessionTraces{1},figOptions};
% % plot_rasters_psth(sessionTraces{1},figOptions);
% 
% sessionTraces{2}=[averages.trialAverages];sessionTraces{2}=sessionTraces{2}(2,1:4:end)';
% sessionTraces{2}=vertcat(sessionTraces{2}{:});
% % sessionTaskRelatedCells{2}=[data.taskRelated];sessionTaskRelatedCells{2}=sessionTaskRelatedCells{2}(2,1:2:end)';
% % sessionTaskRelatedCells{2}=vertcat(sessionTaskRelatedCells{2}{:});
% [~,peakIdx]=max(sessionTraces{2},[],2);
% [~,peakOrder]=sort(peakIdx);
% % sessionTaskRelatedCells{2}=sessionTaskRelatedCells{2}(peakOrder,1);
% sessionTraces{2}=sessionTraces{2}(peakOrder,:);
% 
% figOptions.legends={'Time (s)'; 'Cells'; {'Session 2 (expert) ',...
%         'calcium response aligned to movement initiation'}};
% args{2,1}={sessionTraces{2},figOptions};
% % plot_rasters_psth(sessionTraces{2},figOptions);
% 
% sessionTraces{3}=[averages.trialAverages];sessionTraces{3}=sessionTraces{3}(3,1:4:end)';
% sessionTraces{3}=vertcat(sessionTraces{3}{:});
% % sessionTaskRelatedCells{3}=[data.taskRelated];sessionTaskRelatedCells{3}=sessionTaskRelatedCells{3}(3,1:2:end)';
% % sessionTaskRelatedCells{3}=vertcat(sessionTaskRelatedCells{3}{:});
% [~,peakIdx]=max(sessionTraces{3},[],2);
% [~,peakOrder]=sort(peakIdx);
% % sessionTaskRelatedCells{3}=sessionTaskRelatedCells{3}(peakOrder,1);
% sessionTraces{3}=sessionTraces{3}(peakOrder,:);
% 
% figOptions.legends={'Time (s)'; 'Cells'; {'Session 3 (expert) ',...
%         'calcium response aligned to movement initiation'}};
% args{3,1}={sessionTraces{3},figOptions};
% % plot_rasters_psth(sessionTraces{3},figOptions);
% % plot_rasters_psth(sessionTraces{3}(sessionTaskRelatedCells{3}==1,:),figOptions);
% 
% plot_rasters_psth(args);

% single figure
% averageTracesh=figure; colormap('bone')
% cmapLims=[min(min(vertcat(sessionTraces{:}))), ...
%     max(max(vertcat(sessionTraces{:})))];
% for sessionNum=1:3
%     subplot(1,3,sessionNum)
%     imagesc(sessionTraces{sessionNum},[-2 6]);colorbar;
%     ylabel('Cells','FontWeight','bold','FontSize',12);
%     xlabel('Time (s.)','FontWeight','bold','FontSize',12);
%     currylim=get(gca,'YLim');
%     set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
%     set(gca,'XTickLabel',figOptions.figureStyle{3});
%     % draw alignment bar
%     patch([repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3} repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
%         [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%         [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.3);
%     box off; % axis tight
%     title(['Session ' num2str(sessionNum)]);
% end

%% plot each session for each animal 
% averageTracesh=figure; colormap('bone')
for subjectNum=1:4
    for sessionNum=1:3
        avTraceSubjSessh(subjectNum,sessionNum)=subplot(4,3,sessionNum+3*(subjectNum-1));
        if ~isempty(stats)
            % uncomment section below to create split task-related / non-related figure
%             taskRelatedCells=averages(subjectNum).trialAverages{sessionNum,2}(stats(subjectNum).taskRelated(sessionNum).indices==1,:);
%             [~,peakIdx]=max(taskRelatedCells,[],2);[~,peakOrder]=sort(peakIdx);
%             taskRelatedCells=taskRelatedCells(peakOrder,:);
%             nontaskRelatedCells=averages(subjectNum).trialAverages{sessionNum,2}(stats(subjectNum).taskRelated(sessionNum).indices==0,:);
%             [~,peakIdx]=max(nontaskRelatedCells,[],2);[~,peakOrder]=sort(peakIdx);
%             nontaskRelatedCells=nontaskRelatedCells(peakOrder,:);
%             reorderedRasters=[taskRelatedCells;... %significant
%                                nan(5,70);...        %divider
%                                nontaskRelatedCells]; %non-significant
%             imagesc(reorderedRasters);colorbar;
%         else
            allCells=averages(subjectNum).trialAverages{sessionNum,2};
            [~,peakIdx]=max(allCells,[],2);[~,peakOrder]=sort(peakIdx);
            allCells=allCells(peakOrder,:);
            imagesc(zscore(allCells));colorbar;
        end
        
        ylabel('Cells','FontWeight','bold','FontSize',12);
        currylim=get(gca,'YLim');
        set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
        set(gca,'XTickLabel',figOptions.figureStyle{3});
        % draw alignment bar
        patch([repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3} repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.3);
        box off; % axis tight
        title(['Subject ' num2str(subjectNum) ' Session ' num2str(sessionNum)]);
    end
end

%% isolate the ones with a peak activity around the time of movement
movPeakh=figure;
for subjectNum=1:4
    for sessionNum=1:3
    allCells=averages(subjectNum).trialAverages{sessionNum,2};
    [~,peakIdx]=max(allCells,[],2);[~,peakOrder]=sort(peakIdx);
    allCells=allCells(peakOrder,:);
    % find those with peak around movement time
    movementReponsiveIdx=peakIdx>=16 & peakIdx<=26;
    movementReponsiveIdx=movementReponsiveIdx(peakOrder);
    % frame neurons on previous figure's plots;
    axes(avTraceSubjSessh(subjectNum,sessionNum)); hold on
    patch([repmat(16,1,2) repmat(26,1,2)], ...
        [[find(movementReponsiveIdx, 1 ) find(movementReponsiveIdx, 1, 'last' )]...
        fliplr([find(movementReponsiveIdx, 1 ) find(movementReponsiveIdx, 1, 'last' )])], ...
        [0 0 0 0],'red','EdgeColor','none','FaceAlpha',0.5);
    %plot on 
    figure(movPeakh);
    subplot(4,3,sessionNum+3*(subjectNum-1));
    imagesc(allCells(movementReponsiveIdx,1:50));
    ylabel('Cells','FontWeight','bold','FontSize',12);
    xlabel('Time (s)','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
    set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figOptions.figureStyle{3});
    % draw alignment bars
    patch([repmat(figOptions.alignSpecs{1},1,2) repmat(figOptions.alignSpecs{1}+0.2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    patch([repmat(figOptions.alignSpecs{1},1,2)-5.1 repmat(figOptions.alignSpecs{1}-4.9,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],'red','EdgeColor','none','FaceAlpha',0.5);
    patch([repmat(figOptions.alignSpecs{1},1,2)+4.9 repmat(figOptions.alignSpecs{1}+5.1,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],'red','EdgeColor','none','FaceAlpha',0.5);
    box off; % axis tight
    title(['Subject ' num2str(subjectNum) ' Session' num2str(sessionNum)]);
    end
end

%% same as above, but based on stats
movPeakh=figure;
for subjectNum=1:4
    for sessionNum=1:3
    allCells=averages(subjectNum).trialAverages{sessionNum,2};
    % get movement cells index 
    % use statistatical classification: all positive indices
%     movementReponsiveIdx=stats(subjectNum).taskRelated(sessionNum).indices;

    % use statistatical classification: categories
    movementReponsiveIdx=cellfun(@(class) contains(class,'pre') | contains(class,'during'),...
        stats(subjectNum).taskRelated(sessionNum).classification);

    % use statistatical classification: classify by timing
%     movementReponsiveIdx=arrayfun(@(catCriterion) catCriterion>=-500 && catCriterion<=500,...
%         stats(subjectNum).taskRelated(sessionNum).time); 

    mvtResponseTime=stats(subjectNum).taskRelated(sessionNum).time;
    %sort by peak
    [~,peakIdx]=max(allCells,[],2);[~,peakOrder]=sort(peakIdx);
    allCells=allCells(peakOrder,:);
    movementReponsiveIdx=movementReponsiveIdx(peakOrder);
    mvtResponseTime=mvtResponseTime(peakOrder);
    % frame neurons on previous figure's plots;
    axes(avTraceSubjSessh(subjectNum,sessionNum)); hold on
    plot(figOptions.alignSpecs{1}+mvtResponseTime(movementReponsiveIdx)/100,find(movementReponsiveIdx),'gd')
    %plot on 
    figure(movPeakh);
    subplot(4,3,sessionNum+3*(subjectNum-1));
    imagesc(allCells(movementReponsiveIdx,1:50));
    ylabel('Cells','FontWeight','bold','FontSize',12);
    xlabel('Time (s)','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
    set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figOptions.figureStyle{3});
    % draw alignment bars
    patch([repmat(figOptions.alignSpecs{1},1,2) repmat(figOptions.alignSpecs{1}+0.2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    patch([repmat(figOptions.alignSpecs{1},1,2)-5.1 repmat(figOptions.alignSpecs{1}-4.9,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],'red','EdgeColor','none','FaceAlpha',0.5);
    patch([repmat(figOptions.alignSpecs{1},1,2)+4.9 repmat(figOptions.alignSpecs{1}+5.1,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],'red','EdgeColor','none','FaceAlpha',0.5);
    box off; % axis tight
    title(['Subject ' num2str(subjectNum) ' Session' num2str(sessionNum)]);
    end
end



