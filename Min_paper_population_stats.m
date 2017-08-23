function stats=Min_paper_population_stats(data,params,figOptions)

durationArray=params.durationArray;
preAlignWindow=params.preAlignWindow;
postAlignWindow=params.postAlignWindow;
interval=params.interval;
numSubject=params.numSubject;

stats=struct(data); fldName=fieldnames(stats);
stats=rmfield(stats,fldName(~cellfun(@(x) contains(x,'subject'), fldName)));
[stats.taskRelated, stats.outcomeRelated]=deal(struct('indices',[],'time',[],'peakTime',[],'duration',[],'pVals',[],'classification',{}));

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
        
        for cellNum=1:numCells % stats start here
            %% stats
            %             figure;
            %             allTraces=cell2mat(cellfun(@(mvtimes) caRawTrace(...
            %                 durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
            %                 movementTimes,'UniformOutput', false)');
            allTraces=cell2mat(cellfun(@(mvtimes) caEventTrace(...
                durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
                movementTimes,'UniformOutput', false)');
            %             allTraces=zscore(cell2mat(cellfun(@(mvtimes) caRawTrace(...
            %                 durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
            %                 movementTimes,'UniformOutput', false)'),[],2);
            %             subplot(1,2,1);
            %             imagesc(allTraces);
            %             allTraces=zscore(cell2mat(cellfun(@(mvtimes) caRawTrace(...
            %                 durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
            %                 movementTimes,'UniformOutput', false)'),[],2);
            %             allTraces=zscore(cell2mat(cellfun(@(mvtimes) caEventTrace(...
            %                 durationArray+ceil((mvtimes(1)-preAlignWindow)/interval)-1,cellNum)' ,...
            %                 movementTimes,'UniformOutput', false)'),[],2);
            %             subplot(1,2,2);
            %             imagesc(allTraces);
            
            % Baseline is 1.5-1s before movement.
            %             baselines=mean(allTraces(:,preAlignWindow/interval-15:preAlignWindow/interval-11),2)'; %baseline over 5 frames
            baselines=sum(allTraces(:,1:10),2);
            
            %% Wilcoxon sign rank test (or paired t-test)
            responses=sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5),2);
            [stats(subjectNum).taskRelated(sessionNum).pVals(cellNum,1),...
                stats(subjectNum).taskRelated(sessionNum).indices(cellNum,1)]=...
                signrank(responses,baselines,'Alpha',0.05,'tail','right'); % signrank or ttest
            if stats(subjectNum).taskRelated(sessionNum).indices(cellNum)
                stats(subjectNum).taskRelated(sessionNum).time(cellNum)=preAlignWindow/interval-5+find(sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5))==...
                    max(sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5))),1);
            else
                stats(subjectNum).taskRelated(sessionNum).time(cellNum)=NaN;
            end
            
            %% comparing calcium events frenquency to expected frequency (frequency variant of Wilcoxon test above)
            % However, using ttest (because comparing mean). That's not
            % right. Distributions would need to be log'ed, but then this
            % lead to removing -Inf values, thus only comparing
            % distributions of firing rate when the cell is active
            
            %             baselineMean=sum(reshape(allTraces(:,1:15),[1 size(allTraces,1)*15]))/(size(allTraces,1)*15)*1000;
            %             responses=sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5),2)*100; %it's /10 *1000 for 1s window
            %             [stats(subjectNum).taskRelated(sessionNum).indices(cellNum,1),...
            %                 stats(subjectNum).taskRelated(sessionNum).pVals(cellNum,1)]=ttest(responses,baselineMean,'Alpha',0.01,'tail','right');
            %             if ~isnan(stats(subjectNum).taskRelated(sessionNum).indices(cellNum)) & stats(subjectNum).taskRelated(sessionNum).indices(cellNum)
            %                 stats(subjectNum).taskRelated(sessionNum).time(cellNum)=preAlignWindow/interval-5+find(sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5))==...
            %                     max(sum(allTraces(:,preAlignWindow/interval-4:preAlignWindow/interval+5))),1);
            %             else
            %                 stats(subjectNum).taskRelated(sessionNum).time(cellNum)=NaN;
            %             end
            
            % Kruskal Wallis on bl vs. pre vs. post
            %             preResponses=mean(allTraces(:,preAlignWindow/interval-3:preAlignWindow/interval),2)';
            %             postResponses=mean(allTraces(:,preAlignWindow/interval:preAlignWindow/interval+3),2)';
            %
            %             stats(subjectNum).taskRelated(sessionNum).pVals(cellNum,1)=...
            %                 kruskalwallis([baselines',preResponses',postResponses'],[],'off'); % < 0.01;
            %             stats(subjectNum).taskRelated(sessionNum).indices(cellNum,1)=...
            %             stats(subjectNum).taskRelated(sessionNum).pVals(cellNum,1)<0.05;
            
            %% Bootstrap test
            %             [pVals,direction]=deal(nan(30,1));
            %             for timePoint=-10:19
            %                 timePointReponses=allTraces(:,preAlignWindow/interval+timePoint)';
            %                 [~, ~, pVals(timePoint+11)] = statcond({baselines timePointReponses},...
            %                     'method', 'bootstrap','naccu', 200,'verbose','off');
            %                 direction(timePoint+11)=mean(timePointReponses)-mean(baselines);
            %                 %                                 if pVals(timePoint+11)<=0.01
            %                 %                                     figure; hold on; plot(baselines); plot(timePointReponses);
            %                 %                                 end
            %             end
            %             successiveSignificantVal = regionprops(pVals<=0.01,'Area','PixelIdxList');%find size and index of significant period
            %             stats(subjectNum).taskRelated(sessionNum).pVals(cellNum,:)=pVals;
            %             stats(subjectNum).taskRelated(sessionNum).indices(cellNum)=logical(sum([successiveSignificantVal.Area]>1)) &...
            %                 median(direction(successiveSignificantVal...
            %                 (find([successiveSignificantVal.Area]>1,1)).PixelIdxList))>0; % check that activity is increasing
            %             if stats(subjectNum).taskRelated(sessionNum).indices(cellNum)
            %                 sigDiffTimes=(successiveSignificantVal...
            %                     (find([successiveSignificantVal.Area]>1,1)).PixelIdxList-11)*interval;
            %
            %                 stats(subjectNum).taskRelated(sessionNum).time(cellNum)=sigDiffTimes(1);
            %                 stats(subjectNum).taskRelated(sessionNum).duration(cellNum)=sigDiffTimes(end)-sigDiffTimes(1);
            %
            %                 sigIdx=vertcat(successiveSignificantVal([successiveSignificantVal.Area]>1).PixelIdxList);
            %                 traceExcerpt=allTraces(:,(-10:19)+(preAlignWindow/interval));
            %                 stats(subjectNum).taskRelated(sessionNum).peakTime(cellNum)=...
            %                     (find(mean(traceExcerpt)==max(mean(traceExcerpt)),1)-11)*interval;
            %
            %                 %                 figure; subplot(2,1,1)
            %                 %                 imagesc(traceExcerpt);
            %                 %                 subplot(2,1,2); hold on
            %                 %                 plot(mean(allTraces(:,(-10:19)+(preAlignWindow/interval))))
            %                 %                 plot(sigIdx,ones(numel(sigIdx),1)*max(get(gca,'ylim')),'rd')
            %
            %                 %IMPORTANT: use either sigDiffTimes(1) or peakTime to categorize neurons
            %                 % peakTime is the time of maximum value within the 3s task period (-1 to +2s),
            %                 % NOT within the significant time points. Which can lead to a peakTime
            %                 % earlier than the sigDiffTimes(1) :(
            %                 catCriterion= sigDiffTimes(1); % stats(subjectNum).taskRelated(sessionNum).peakTime(cellNum); %sigDiffTimes(1);
            %                 if catCriterion>=-1000 && catCriterion<-500
            %                     stats(subjectNum).taskRelated(sessionNum).classification{cellNum,1}='early';
            %                 elseif catCriterion>=-500 && catCriterion<0
            %                     stats(subjectNum).taskRelated(sessionNum).classification{cellNum,1}='pre';
            %                 elseif catCriterion>=0 && catCriterion<500
            %                     stats(subjectNum).taskRelated(sessionNum).classification{cellNum,1}='during';
            %                 elseif catCriterion>=500 && catCriterion<2000
            %                     stats(subjectNum).taskRelated(sessionNum).classification{cellNum,1}='late';
            %                 end
            %
            %                 %% example stat plot (similar to Komiyama 2010 Fig S10)
            % %                                 extra figure options
            % %                                 figOptions.figureStyle={'hot';1:10:size(allTraces,2)+1;-preAlignWindow/1000:postAlignWindow/1000};
            % %                                 figOptions.legends={'Time (s)'; 'Trials'; {['Neuron # ' num2str(cellNum)],...
            % %                                     'calcium response aligned to movement initiation'}};
            % %                                 figure; colormap(figOptions.figureStyle{1})
            % %                                 plot rasters
            % %                                 subplot(2,2,[1,3])
            % %                                 imagesc(allTraces); colorbar;
            % %                                 ylabel('Trials','FontWeight','bold','FontSize',12);
            % %                                 currylim=get(gca,'YLim');
            % %                                 set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
            % %                                 set(gca,'XTickLabel',figOptions.figureStyle{3});
            % %                                 draw alignment bar
            % %                                 patch([repmat(figOptions.alignSpecs{1},1,2) repmat(figOptions.alignSpecs{1}+0.5,1,2)], ...
            % %                                     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            % %                                     [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
            % %                                 box off;
            % %
            % %                                 plot mean trace
            % %                                 subplot(2,2,2); hold on
            % %                                 sem=std(allTraces)/ sqrt(size(allTraces,1));
            % %                                 patch([1:length(sem),fliplr(1:length(sem))],[mean(allTraces)-sem,...
            % %                                     fliplr(mean(allTraces)+sem)],[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.3);
            % %                                 plot(mean(allTraces),'LineWidth',1.8);%'Color',cmap(rastnum,:)
            % %                                 currylim=get(gca,'YLim');
            % %                                 set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
            % %                                 set(gca,'XTickLabel',figOptions.figureStyle{3});
            % %                                 draw alignment bar
            % %                                 patch([repmat(figOptions.alignSpecs{1},1,2) repmat(figOptions.alignSpecs{1}+0.5,1,2)], ...
            % %                                     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            % %                                     [0 0 0 0],'k','EdgeColor','none','FaceAlpha',0.5);
            % %                                 title('Mean calcium trace')
            % %                                 box off; axis tight
            % %
            % %                                 plot p values
            % %                                 subplot(2,2,4); hold on
            % %                                 plot(pVals,'LineWidth',1.8);
            % %                                 plot(vertcat(successiveSignificantVal.PixelIdxList),...
            % %                                     ones(length(vertcat(successiveSignificantVal.PixelIdxList)),1)*0.2,'k*');
            % %                                 set(gca,'ylim',[0 0.5]);
            % %                                 currylim=get(gca,'YLim');
            % %                                 set(gca,'XTick',[1,30],'TickDir','out');
            % %                                 set(gca,'XTickLabel',round(linspace(-1000,1900,2)));
            % %                                 draw alignment bar
            % %                                 patch([repmat(10,1,2) repmat(10.5,1,2)], ...
            % %                                     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            % %                                     [0 0 0 0],'k','EdgeColor','none','FaceAlpha',0.5);
            % %                                 title('p values')
            % %                                 box off;
            %
            %             else
            %                 stats(subjectNum).taskRelated(sessionNum).time(cellNum)=NaN;
            %                 stats(subjectNum).taskRelated(sessionNum).duration(cellNum)=NaN;
            %                 stats(subjectNum).taskRelated(sessionNum).classification{cellNum,1}='nonresponsive';
            %             end
            
            %% kruskalwallis test
            % this is the version when epochs are 5 frames + movement
            % duration. But it's too complicated to deal with uneven bin
            % numbers
            %             binedMovementEpochs =  cellfun(@(trialMvTime) accumarray(floor(linspace(1,ceil(numel(trialMvTime)/2),...
            %                 numel(caEvents(cellNum,trialMvTime))))',caEvents(cellNum,trialMvTime),[],@mean), ...
            %                 cellfun(@(trialMvTime) preAlignWindow/100-5:(preAlignWindow/100+ceil(diff(trialMvTime)/100)),...
            %                 movementTimes,'UniformOutput', false),'UniformOutput', false);
            % simplified version
            %             trialMvArray=cellfun(@(trialArray)(trialArray.Time(preAlignWindow/100-5:preAlignWindow/100+2)-1)/100,...
            %                 {data(subjectNum).ciData(trialIdx).rawTraceEpochs},'UniformOutput', false);
            %             binedMovementEpochs = cellfun(@(trialMvTime) rot90(accumarray(floor(linspace(1,ceil(numel(trialMvTime)/2),...
            %                 numel(caRawTrace.Data(trialMvTime,cellNum))))',caRawTrace.Data(trialMvTime,cellNum),[],@mean)), ...
            %                 trialMvArray,'UniformOutput', false);
            %             binedMovementEpochs = cellfun(@(trialMvTime) (accumarray(sort(repmat(1:4,1,2))',caRawTrace(...
            %                 (0:7)+ceil((trialMvTime(1)-500)/interval),cellNum),[],@mean))', ... %get 8 bins, starting 500ms before movement
            %                 movementTimes,'UniformOutput', false);
            
            % working but stringent response test
            %             binedMovementEpochs = cellfun(@(trialMvTime) (caRawTrace(...
            %                 (0:7)+ceil((trialMvTime(1)-500)/interval),cellNum))', ... %get 8 bins, starting 500ms before movement
            %                 movementTimes,'UniformOutput', false);
            %
            %             stats(subjectNum).taskRelated{sessionNum}(cellNum) = kruskalwallis(vertcat(binedMovementEpochs{:}),[],'off'); % < 0.01;
            
            
            %% now split and compare trials by outcome (success, fail)
            %             compare %m:missed vs f:failed (may add droped)
            %                         pVals=nan(10,1);
            %                         for timePoint=-5:4
            %                             successTrialsReponses=allTraces(outcomeIdx.success,preAlignWindow/interval+timePoint)';
            %                             failedTrialsReponses=allTraces(outcomeIdx.failed,preAlignWindow/interval+timePoint)';
            %                             [~, ~, pVals(timePoint+6)] = statcond({successTrialsReponses failedTrialsReponses},...
            %                                 'method', 'bootstrap','naccu', 2000,'verbose','off');
            %                             %                 if pVals(timePoint+1)<=0.01
            %                             %                     figure; hold on; plot(baselines); plot(timePointReponses);
            %                             %                 end
            %                         end
            %                         successiveSignificantVal = regionprops(pVals<=0.01,'Area','PixelIdxList');%find size and index of significant period
            %                         stats(subjectNum).outcomeRelated(sessionNum).pVals(cellNum,:)=pVals;
            %                         stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum)=logical(sum([successiveSignificantVal.Area]>1));
            %                         if stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum)
            %                             successiveSignificantVal=successiveSignificantVal([successiveSignificantVal.Area]>1);
            %                             sigDiffTimes=(successiveSignificantVal(1).PixelIdxList-6)*interval; %keeping first "significant" time window
            %                             stats(subjectNum).outcomeRelated(sessionNum).time(cellNum)=sigDiffTimes(1);
            %                             stats(subjectNum).outcomeRelated(sessionNum).duration(cellNum)=sigDiffTimes(end)-sigDiffTimes(1);
            %                             sigIdx=preAlignWindow/interval+reshape(vertcat(successiveSignificantVal.PixelIdxList),...
            %                                 [numel(vertcat(successiveSignificantVal.PixelIdxList)),1])-6; %Here all significant period taken into account
            %                             if max(mean(allTraces(outcomeIdx.success,sigIdx)))>...
            %                                     max(mean(allTraces(outcomeIdx.failed,sigIdx)))
            %                                 stats(subjectNum).outcomeRelated(sessionNum).classification{cellNum,1}='success';
            %                             else
            %                                 stats(subjectNum).outcomeRelated(sessionNum).classification{cellNum,1}='failure';
            %                             end
            %                         else
            %                             stats(subjectNum).outcomeRelated(sessionNum).time(cellNum)=NaN;
            %                             stats(subjectNum).outcomeRelated(sessionNum).duration(cellNum)=NaN;
            %                             stats(subjectNum).outcomeRelated(sessionNum).classification{cellNum,1}='non_outcome_related';
            %                         end
            
            %% test movement response vs outcome
            %test success vs. failure, success vs. miss, success vs. no seed
            respVsOutcome=cell2table(trialOutcomes','VariableNames',{'trialOutcomes'});
            respVsOutcome.responses=sum(allTraces(:,16:25),2)>0;
            % success vs. failure
            respVsOutcome=respVsOutcome(contains(respVsOutcome.trialOutcomes,'f')|...
                contains(respVsOutcome.trialOutcomes,'s'),:);
            contingencyTbl = crosstab(respVsOutcome.trialOutcomes,respVsOutcome.responses);
            try
                [stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,1),...
                    stats(subjectNum).outcomeRelated(sessionNum).pVals(cellNum,1)] = fishertest(contingencyTbl);
            catch
                stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,1) = 0;
            end
            % success vs. miss
            respVsOutcome=respVsOutcome(contains(respVsOutcome.trialOutcomes,'m')|...
                contains(respVsOutcome.trialOutcomes,'s'),:);
            contingencyTbl = crosstab(respVsOutcome.trialOutcomes,respVsOutcome.responses);
            try
                [stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,2),...
                    stats(subjectNum).outcomeRelated(sessionNum).pVals(cellNum,2)] = fishertest(contingencyTbl);
            catch
                stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,2) = 0;
            end
            % success vs. no seed
            respVsOutcome=respVsOutcome(contains(respVsOutcome.trialOutcomes,'n')|...
                contains(respVsOutcome.trialOutcomes,'s'),:);
            contingencyTbl = crosstab(respVsOutcome.trialOutcomes,respVsOutcome.responses);
            try
                [stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,3),...
                    stats(subjectNum).outcomeRelated(sessionNum).pVals(cellNum,3)] = fishertest(contingencyTbl);
            catch
                stats(subjectNum).outcomeRelated(sessionNum).indices(cellNum,3) = 0;
            end
        end
        
    end
end





