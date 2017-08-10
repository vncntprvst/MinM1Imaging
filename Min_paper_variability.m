function Min_paper_variability(data,stats,params)

for subjectNum=1:4
    % here using the statistical category - can be changed 
    %         movementReponsiveIdx=cellfun(@(class) contains(class,'pre') | contains(class,'during'),...
    %             stats(subjectNum).taskRelated(sessionNum).classification);
    movementReponsiveIdx=arrayfun(@(catCriterion) catCriterion>=-500 && catCriterion<=500,...
        stats(subjectNum).taskRelated(1).time);
    
    % forward destiny: what Naive Session movement cell become
    cellIDs=data(subjectNum).cellIDs(1).cellIndex;
    [~,earlySessionIdx]=sort(cellIDs.Early);
    cellIDs=cellIDs(earlySessionIdx,:);cellIDs=cellIDs(cellIDs.Early~=0,:);
    
    movementReponsiveCells=find(movementReponsiveIdx);
    mvtCellsIDs{1}=movementReponsiveCells(ismember(movementReponsiveCells,cellIDs.Early));
    
    mvtCellsIDs{2}=cellIDs.Expert1(ismember(cellIDs.Early,mvtCellsIDs{1}));
    mvtCellsIDs{3}=cellIDs.Expert2(ismember(cellIDs.Early,mvtCellsIDs{1}));
    
    % find bad / no spike / too short trials
    badTrials=cellfun(@(spikeCell) size(spikeCell,2)==1,...
        {data(subjectNum).ciData.spikes})'; %likely interrupted trial at the end of a session
    sampleLength=cell2mat(cellfun(@(trialResponse) get(trialResponse,'Length'),...
        {data(subjectNum).ciData.rawTraceEpochs}','UniformOutput', false));
    shortTrials=sampleLength<max(params.durationArray); %to exclude trials that are too short
    sessionDays=unique([data(subjectNum).behavData.session]);
    
    % spike timing synchrony
    for sessionNum=1:3
        %indices
        sessionIdx=[data(subjectNum).behavData.session]'==sessionDays(sessionNum);
        trialIdx=~(badTrials | shortTrials | ~sessionIdx);
        
        % get traces
        caRawTrace=data(subjectNum).rawTraces{sessionNum,1}.Data;
        caEventTrace=data(subjectNum).calciumEvents{sessionNum};
        movementTimes={data(subjectNum).behavData(trialIdx).movementTime};
        
        % get t trials cell, with 70*n neurons in each cell 
        mvtNeuronsTraces=cellfun(@(mvtimes) caRawTrace(params.durationArray+ceil((mvtimes(1)-...
            params.preAlignWindow)/params.interval)-1,mvtCellsIDs{sessionNum}),...
            movementTimes,'UniformOutput', false);
        
        % need to normalize data: 
        %         For analysing the longitudinal dynamics of the fractions of movement-related
        %         neurons over sessions (Fig. 2b), the fraction of movement-related neurons in each
        %         session in each animal was normalized to the highest and lowest values of the animal.
        
%         correlations{sessionNum}=cellfun(@(trialData) xcorr(trialData'), mvtNeuronsTraces,'UniformOutput', false);
%         correlations{sessionNum}=cellfun(@(trialData) corr(trialData), mvtNeuronsTraces,'UniformOutput', false);
%         correlations{sessionNum}=cellfun(@(trialData) corrcoef(trialData), mvtNeuronsTraces,'UniformOutput', false);
        
        correlations{sessionNum}=cellfun(@(trialData) corr(trialData'), mvtNeuronsTraces,'UniformOutput', false);

        figure; imagesc(correlations{1, 1}{1, 70})
        
        foo=mean(cat(3,correlations{sessionNum}{:}),3);  
        figure; imagesc(foo);
        
        % get t trials cell, with 20 time points * n neurons in each cell 
        mvtNeuronsTraces=cellfun(@(mvtimes) caRawTrace((1:20)+ceil((mvtimes(1)-...
            1000)/params.interval)-1,mvtCellsIDs{sessionNum}),...
            movementTimes,'UniformOutput', false);

        correlations{sessionNum}=cellfun(@(trialData) corr(trialData'), mvtNeuronsTraces,'UniformOutput', false);

        figure; imagesc(correlations{1, 1}{1, 30})
        
        foo=mean(cat(3,correlations{sessionNum}{:}),3);  
        figure; imagesc(foo)
        

    end
end

