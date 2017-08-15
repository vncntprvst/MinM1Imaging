function [meanTimeCorrelation,meanTimeCovariance]=Min_paper_variability(data,stats,params,figOptions)

for subjectNum=1:4
    %% forward destiny: what Naive Session movement cell become
    %  here using the statistical category - can be changed
    %     %         movementReponsiveIdx=cellfun(@(class) contains(class,'pre') | contains(class,'during'),...
    %     %             stats(subjectNum).taskRelated(sessionNum).classification);
    %     movementReponsiveIdx=arrayfun(@(catCriterion) catCriterion>=-500 && catCriterion<=500,...
    %         stats(subjectNum).taskRelated(1).time);
    %
    %
    %     cellIDs=data(subjectNum).cellIDs(1).cellIndex;
    %     [~,earlySessionIdx]=sort(cellIDs.Early);
    %     cellIDs=cellIDs(earlySessionIdx,:);cellIDs=cellIDs(cellIDs.Early~=0,:);
    %
    %     movementReponsiveCells=find(movementReponsiveIdx);
    %     mvtCellsIDs{1}=movementReponsiveCells(ismember(movementReponsiveCells,cellIDs.Early));
    %     mvtCellsIDs{2}=cellIDs.Expert1(ismember(cellIDs.Early,mvtCellsIDs{1}));
    %     mvtCellsIDs{2}=mvtCellsIDs{2}(mvtCellsIDs{2}~=0);
    %     mvtCellsIDs{3}=cellIDs.Expert2(ismember(cellIDs.Early,mvtCellsIDs{1}));
    %     mvtCellsIDs{3}=mvtCellsIDs{3}(mvtCellsIDs{3}~=0);
    
    %% find bad / no spike / too short trials
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
        
        %         cellIDs= mvtCellsIDs{sessionNum};  %if doing destiny
        cellIDs= arrayfun(@(catCriterion) catCriterion>=-500 && catCriterion<=500,...
            stats(subjectNum).taskRelated(sessionNum).time);
        
        % get t trials cell, with 70*n neurons in each cell
        mvtNeuronsTraces=cellfun(@(mvtimes) caRawTrace(params.durationArray+ceil((mvtimes(1)-...
            params.preAlignWindow)/params.interval)-1,cellIDs),...
            movementTimes,'UniformOutput', false);
        % get equivalent traces for non-movement cells
        numMvtCells=sum(cellIDs);
        cellIDs=find(~cellIDs); 
         % keep only an equivalent number of cells to mvt cells, randomly selected
         %cellIDs=cellIDs(randperm(numel(cellIDs),numMvtCells));
        nMvtNeuronsTraces=cellfun(@(mvtimes) caRawTrace(params.durationArray+ceil((mvtimes(1)-...
            params.preAlignWindow)/params.interval)-1,cellIDs),...
            movementTimes,'UniformOutput', false);
        
        % need to normalize data:
        %         For analysing the longitudinal dynamics of the fractions of movement-related
        %         neurons over sessions (Fig. 2b), the fraction of movement-related neurons in each
        %         session in each animal was normalized to the highest and lowest values of the animal.
        
        %         correlations{sessionNum}=cellfun(@(trialData) xcorr(trialData'), mvtNeuronsTraces,'UniformOutput', false);
        %         correlations{sessionNum}=cellfun(@(trialData) corr(trialData), mvtNeuronsTraces,'UniformOutput', false);
        %         correlations{sessionNum}=cellfun(@(trialData) corrcoef(trialData), mvtNeuronsTraces,'UniformOutput', false);
        
        mvtNeuronTimeCorrelations{subjectNum,sessionNum}=cellfun(@(trialData) corr(trialData'), mvtNeuronsTraces,'UniformOutput', false);
        mvtNeuronTimeCovariance{subjectNum,sessionNum}=cellfun(@(trialData) cov(trialData'), mvtNeuronsTraces,'UniformOutput', false);
        nMvtNeuronTimeCovariance{subjectNum,sessionNum}=cellfun(@(trialData) cov(trialData'), nMvtNeuronsTraces,'UniformOutput', false);

        %% get mean correlation and covariance matrix
        %         figure; imagesc(mvtNeuronTimeCorrelations{1, 1}{1, 26})
        
        meanTimeCorrelation{subjectNum,sessionNum}=mean(cat(3,mvtNeuronTimeCorrelations{subjectNum,sessionNum}{:}),3);
        %         figure; imagesc(meanTimeCorrelation{sessionNum});
        meanTimeCovariance{subjectNum,sessionNum}=mean(cat(3,mvtNeuronTimeCovariance{subjectNum,sessionNum}{:}),3);
        %         figure; imagesc(meanTimeCovariance{subjectNum,sessionNum});
        meanNMNTimeCovariance{subjectNum,sessionNum}=mean(cat(3,nMvtNeuronTimeCovariance{subjectNum,sessionNum}{:}),3);
               
        %% tried to separate correlation matrices by profile. Didn't work much
        %         corrMatStdProfile=cell2mat(cellfun(@(corrMat) nanstd(corrMat),...
        %             mvtNeuronTimeCorrelations{1, 1},'UniformOutput', false)');
        %
        %         Y = pdist(corrMatStdProfile);
        %         Z = linkage(Y,'average');
        %         %         dendrogram(Z);
        %         c = cophenet(Z,Y)
        %         dendrogram(Z)
        %         c = cluster(Z,'maxclust',3);
        %
        %         figure;
        %         for clusNum=1:max(c)
        %             subplot(3,1,clusNum);
        %             plot(mean(corrMatStdProfile(c==clusNum,:)));
        %         end
        %
        %         figure; clusNum=find(c==1);
        %         for corrMat=1:sum(c==1)
        %             subplot(10,10,corrMat)
        %             imagesc(mvtNeuronTimeCorrelations{1, 1}{1, clusNum(corrMat)});
        %         end
        %         figure; clusNum=find(c==2);
        %         for corrMat=1:sum(c==2)
        %             subplot(10,10,corrMat)
        %             imagesc(mvtNeuronTimeCorrelations{1, 1}{1, clusNum(corrMat)});
        %         end
        %         figure; clusNum=find(c==3);
        %         for corrMat=1:sum(c==3)
        %             subplot(10,10,corrMat)
        %             imagesc(mvtNeuronTimeCorrelations{1, 1}{1, clusNum(corrMat)});
        %         end
        
        %% focus on -1s +1s around movement
        %         % get t trials cell, with 20 time points * n neurons in each cell
        %         mvtNeuronsTraces=cellfun(@(mvtimes) caRawTrace((1:20)+ceil((mvtimes(1)-...
        %             1000)/params.interval)-1,mvtCellsIDs{sessionNum}),...
        %             movementTimes,'UniformOutput', false);
        %
        %         correlations{sessionNum}=cellfun(@(trialData) corr(trialData'), mvtNeuronsTraces,'UniformOutput', false);
        %
        %         figure; imagesc(correlations{1, 1}{1, 30})
        %
        %         foo=mean(cat(3,correlations{sessionNum}{:}),3);
        %         figure; imagesc(foo)
        %
        
    end
end
%% Mean covariance matrix, per session. Pool subjects together
% Movement cells
mvtCellCov{1}=mean(cat(3,meanTimeCovariance{:,1}),3);
mvtCellCov{2}=mean(cat(3,meanTimeCovariance{:,2}),3);
mvtCellCov{3}=mean(cat(3,meanTimeCovariance{:,3}),3);
figure; figOptions.figureStyle={'parula';1:10:size(params.durationArray,2)+1;...
    -params.preAlignWindow/1000:params.postAlignWindow/1000};
cmap=colormap;
for covplot=1:3
    covVals=mvtCellCov{covplot}; covVals=covVals(1:40,1:40);
    subplot(2,3,covplot);
    imagesc(covVals);
    ylabel('Time (s.)','FontWeight','bold','FontSize',12);
    xlabel('Time (s.)','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
%     set(gca,'XTick',figOptions.figureStyle{2},'YTick',figOptions.figureStyle{2},'TickDir','out');
%     set(gca,'XTickLabel',figOptions.figureStyle{3},'YTickLabel',figOptions.figureStyle{3});  
    % draw alignment bar
    patch([repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3} repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    patch([[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3}...
        repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
        figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    box off; % axis tight
    title([' Session ' num2str(covplot)]);
    
    % plot average and sem
    subplot(2,3,covplot+3);hold on
    plot(mean(covVals));
    ylabel('Mean Covariance','FontWeight','bold','FontSize',12);
    xlabel('Time (s.)','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
%     set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
%     set(gca,'XTickLabel',figOptions.figureStyle{3}); 
    % draw sem
        patch([1:size(covVals,1) flip(1:size(covVals,1))], ...
        [mean(covVals)+(std(covVals)/sqrt(size(covVals,1)))...
        fliplr(mean(covVals)-(std(covVals)/sqrt(size(covVals,1))))], ...
        cmap(1,:),'EdgeColor','none','FaceAlpha',0.3);
    % draw alignment bar
    patch([repmat(figOptions.alignSpecs{1},1,2)-0.2 ...
        repmat(figOptions.alignSpecs{1}+0.2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        'k','EdgeColor','none','FaceAlpha',0.5);
    box off; axis tight
end

% non-movement cells
nMvtCellCov{1}=mean(cat(3,meanNMNTimeCovariance{:,1}),3);
nMvtCellCov{2}=mean(cat(3,meanNMNTimeCovariance{:,2}),3);
nMvtCellCov{3}=mean(cat(3,meanNMNTimeCovariance{:,3}),3);
figure; figOptions.figureStyle={'parula';1:10:size(params.durationArray,2)+1;...
    -params.preAlignWindow/1000:params.postAlignWindow/1000};
for covplot=1:3
    covVals=nMvtCellCov{covplot};covVals=covVals(1:40,1:40);
    subplot(2,3,covplot);
    imagesc(covVals);
    ylabel('Time','FontWeight','bold','FontSize',12);
    xlabel('Time','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
%     set(gca,'XTick',figOptions.figureStyle{2},'YTick',figOptions.figureStyle{2},'TickDir','out');
%     set(gca,'XTickLabel',figOptions.figureStyle{3},'YTickLabel',figOptions.figureStyle{3});  
    % draw alignment bar
    patch([repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3} repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    patch([[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [repmat(figOptions.alignSpecs{1},1,2)-figOptions.alignSpecs{3}...
        repmat(figOptions.alignSpecs{1}+figOptions.alignSpecs{3},1,2)], ...
        [0 0 0 0],figOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.5);
    box off; % axis tight
    title([' Session ' num2str(covplot)]);
        
    % plot average and sem
    subplot(2,3,covplot+3);hold on
    plot(mean(covVals));
    ylabel('Mean Covariance','FontWeight','bold','FontSize',12);
    xlabel('Time (s.)','FontWeight','bold','FontSize',12);
    currylim=get(gca,'YLim');
%     set(gca,'XTick',figOptions.figureStyle{2},'TickDir','out');
%     set(gca,'XTickLabel',figOptions.figureStyle{3}); 
    % draw sem
        patch([1:size(covVals,1) flip(1:size(covVals,1))], ...
        [mean(covVals)+(std(covVals)/sqrt(size(covVals,1)))...
        fliplr(mean(covVals)-(std(covVals)/sqrt(size(covVals,1))))], ...
        cmap(1,:),'EdgeColor','none','FaceAlpha',0.3);
    % draw alignment bar
    patch([repmat(figOptions.alignSpecs{1},1,2)-0.2 ...
        repmat(figOptions.alignSpecs{1}+0.2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        'k','EdgeColor','none','FaceAlpha',0.5);
    box off; axis tight
end
