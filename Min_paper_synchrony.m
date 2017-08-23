function syncVals=Min_paper_synchrony(data,stats,params,figOptions)
% spike timing synchrony
for figNum=1:10
    figh{figNum}=figure;
end
% get traces from all movement related cells, for each trial
for subjectNum=1:4
    % get traces
    caRawTraces=cellfun(@(traces,mvtcellsidx) traces.Data(:,mvtcellsidx) ,...
        data(subjectNum).rawTraces(:,1),{stats(subjectNum).taskRelated.indices}','UniformOutput',false);
    
    % find bad / no spike / too short trials
    badTrials=cellfun(@(spikeCell) size(spikeCell,2)==1,...
        {data(subjectNum).ciData.spikes})'; %likely interrupted trial at the end of a session
    sampleLength=cell2mat(cellfun(@(trialResponse) get(trialResponse,'Length'),...
        {data(subjectNum).ciData.rawTraceEpochs}','UniformOutput', false));
    shortTrials=sampleLength<max(params.durationArray); %to exclude trials that are too short
    sessionDays=unique([data(subjectNum).behavData.session]);
    for sessionNum=1:numel(sessionDays)
        sessionIdx=[data(subjectNum).behavData.session]'==sessionDays(sessionNum);
        trialIdx=~(badTrials | shortTrials | ~sessionIdx);
        movementTimes={data(subjectNum).behavData(trialIdx).movementTime};
        
        %get trial traces
        allTraces=cellfun(@(mvtimes) caRawTraces{sessionNum}(...
            params.durationArray+ceil((mvtimes(1)-params.preAlignWindow)/params.interval)-1,:)' ,...
            movementTimes,'UniformOutput', false)';
        
        %         traceCorr=cell2mat(cellfun(@(traces) mean(corr(traces)), allTraces,'UniformOutput', false));
        %         corr(mean(allTraces{1}(:,1:40))',allTraces{1}(1,1:40)')
        
        %         foo=(corr(mean(allTraces{1}(:,1:40))',allTraces{1}(1,1:40)'))
        %         figure; imagesc(foo)
        
        
        %% simplest way: correlations
        %         window=1;
        %         timeCoefs=nan(size(allTraces,1),30);
        %         for steps=11:40
        %             tracesCorr=cellfun(@(trace) corrcoef(trace(:,steps-window:steps+window-1)'),allTraces,'UniformOutput', false);
        %             timeCoefs(:,steps-10)=cell2mat(cellfun(@(ccoefs) mean(abs(ccoefs(logical(tril(ccoefs,-1))))),tracesCorr,'UniformOutput', false));
        %         end
        %         subplot(4,3,sessionNum+3*(subjectNum-1));
        %         imagesc(timeCoefs)
        
        %% semblance
        % mean population vs each cells
        timeSemblance=nan(size(allTraces,1),40,10);
        for trialNum=1:size(allTraces,1)
            traces=allTraces{trialNum};
            trialTimeSemblance=nan(size(traces,1),40,10);
            for cellNum=1:size(traces,1)
                trialTimeSemblance(cellNum,:,:)=semblance(mean(traces(:,1:40)),traces(cellNum,1:40),10)';
            end
            timeSemblance(trialNum,:,:)=mean(trialTimeSemblance,1);
        end
        for figNum=1:10
            figure(figh{figNum});
            subplot(4,3,sessionNum+3*(subjectNum-1));
            imagesc(timeSemblance(:,:,figNum))
        end
        
        %% coherence with wavelets
        % calculate coherence on wavelets (instead of fourrier)
        
        % wavelets x convolved signal
        % coherence on these values
        
        % or use Singular Spectrum Decomposition
        
    end
end







