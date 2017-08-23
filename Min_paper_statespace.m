function Min_paper_statespace(data,stats,params,figOptions)

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
        
        %         [mapped_data, mapping] = compute_mapping(allTraces{1}, 'tSNE'); %, # of dimensions, parameters)
        
        figure;
        cmap=colormap('parula');
        
        rng('default');
        lowDimData = tsne(allTraces{1},'Algorithm','exact','Distance','cosine');
        distLinkage = linkage(lowDimData,'single','cosine');
        cellClusters = cluster(distLinkage,'maxclust',2);
        subplot(2,2,1)
        [~,foo]=sort(lowDimData(:,1));
        gscatter(lowDimData(foo,1),lowDimData(foo,2),1:54,cmap(1:54,:)); %cellClusters
        
        title('Cosine')
        
        rng('default');
        lowDimData = tsne(allTraces{1},'Algorithm','exact','Distance','cityblock'); %'Standardize',true,'Perplexity',20);
        distLinkage = linkage(lowDimData,'single','cityblock');
        cellClusters = cluster(distLinkage,'maxclust',3);
        subplot(2,2,2)
        gscatter(lowDimData(foo,1),lowDimData(foo,2),1:54,cmap(1:54,:)); %cellClusters
        
        title('City block')
        
        rng('default');
        lowDimData = tsne(allTraces{1},'Algorithm','exact','Distance','chebychev');
        distLinkage = linkage(lowDimData,'single','chebychev');
        cellClusters = cluster(distLinkage,'maxclust',3);
        subplot(2,2,3)
        gscatter(lowDimData(foo,1),lowDimData(foo,2),1:54,cmap(1:54,:)); %cellClusters
        title('Chebychev')
        
        rng('default');
        lowDimData = tsne(allTraces{1},'Algorithm','exact','Distance','euclidean');
        distLinkage = linkage(lowDimData,'ward','euclidean');
        cellClusters = cluster(distLinkage,'maxclust',3);
        subplot(2,2,4)
        gscatter(lowDimData(foo,1),lowDimData(foo,2),1:54,cmap(1:54,:)); %cellClusters
        title('Euclidean')
        
        
    end
    
end
