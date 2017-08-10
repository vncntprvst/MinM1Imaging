function [ciData,rawTraces,calciumEvents,cellIDs,behavData]=Min_paper_getCIdata(exportFolder,trialWindow,samplingRate,subject)

% (earlier step suggestion:
% exclude neurons with ratio of positive to negative transients below ten.
% (see Komiyama et al 2010)

% raw data folder
rawDataFolder='D:\Data\Raw\CI\';
cd([rawDataFolder subject filesep 'ana']);
% cd('ana20161023')

%define parameters
if nargin<2
    preAlignWindow=2000; % 2 seconds before rising phase
    postAlignWindow=5000; % 5 seconds after peak
    samplingRate=10; %10Hz
else
    preAlignWindow=trialWindow(1);
    postAlignWindow=trialWindow(2);
end
interval=1000/samplingRate;

dataDirListing=dir;
%keep data files
dataDirListing=dataDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'.mat'),...
    {dataDirListing.name},'UniformOutput',false)));
% subject=cell2mat(unique(cellfun(@(x) x(1:3),{dataDirListing.name},'UniformOutput',false)));
days=cellfun(@(x) x(numel(subject)+1:numel(subject)+3),{dataDirListing.name},'UniformOutput',false);
daysIdx=contains(days,'d');
days=unique(cellfun(@(x,y) str2double(x(~isletter(regexprep(y(numel(subject)+1:numel(subject)+3),'\.','')))),...
    days(daysIdx),{dataDirListing(daysIdx).name},'UniformOutput',true));

% get behavior data
MinBehavData=load([subject 'beh.mat']);
% structure definition
behavData=struct('subject',[],'session',[],'trialNum',[],...
    'movementTime',[],'outcome',[]);
ciData=struct('trialWindowIndex',[],'alignTime',[],'rawTraceEpochs',[],'spikes',[]);
cellIDs=struct('cellIndex',[],'ROIs',[]);
rawTraces=cell(numel(days),2);calciumEvents=cell(numel(days),1);
trialId=1;
for session=1:numel(days)
    try
        MinCIData=load([subject 'd' num2str(days(session)) '.mat']);
        % keep raw signals
        rawTraces{session,1}=[MinCIData.([subject 'd' num2str(days(session))]){2:end, 4}];
        rawTraces{session,1}=timeseries(rawTraces{session},1:interval:(interval)*numel(MinCIData.([subject 'd' num2str(days(session))]){session+1, 4}),...
            'Name','Calcium signal');
        % event detection: 
            % three times the median absolute deviation (MAD))
            % then non-negative deconvolution method - Huber 2012
        g = 0.92; lambda=2.4;       
        rawTraces{session,2}=rawTraces{session,1};
        traces=rawTraces{session,2}.Data;
        [c_oasis_ar1, s_oasis_ar1] = deconvolveCa(traces, 'ar1', g, 'foopsi', 'lambda', lambda);
        rawTraces{session,2}.Data=reshape(c_oasis_ar1,[],size(traces,2));
        calciumEvents{session}=reshape(s_oasis_ar1,[],size(traces,2));   

        ROIs=[MinCIData.([subject 'd' num2str(days(session))]){2:end, 2}];
        % keep ROIs and cell index
        cellIDs(session).ROIs=[ROIs(1:2:end)' ROIs(2:2:end)'];
        cellIDs(session).cellIndex=importdata([subject 'cellIDs.mat']);
        timeInfo=get(rawTraces{session},'Timeinfo');
        allSpikes={MinCIData.([subject 'd' num2str(days(session))]){2:end, 7}};
        allSpikes = cellfun(@(spikeTimes) spikeTimes*interval, allSpikes, 'UniformOutput', false);
        numTrials=numel(MinBehavData.([subject 'beh']){session+1, 2});
        for trialNum=trialId:trialId+numTrials-1
            behavData(trialNum).subject=subject;
            behavData(trialNum).session=str2double(num2str(days(session)));
            behavData(trialNum).trialNum=trialNum-trialId+1;
            behavData(trialNum).movementTime=round(MinBehavData.([subject 'beh'])...
                {session+1, 4}(behavData(trialNum).trialNum,:)*(interval)); %in ms
            behavData(trialNum).outcome=MinBehavData.([subject 'beh'])...
                {session+1, 2}(behavData(trialNum).trialNum,:);
            
            % import the fluorescence trace
            %         ciData(trialNum).session=str2double(num2str(days(session)));
            %         ciData(trialNum).trialNum=trialNum-trialId+1;
            ciData(trialNum).trialWindowIndex=[max([round(behavData(trialNum).movementTime(1))-preAlignWindow 1]) ...
                min([round(behavData(trialNum).movementTime(2))+postAlignWindow timeInfo.End])];
            ciData(trialNum).alignTime=[max([round(behavData(trialNum).movementTime(1))-...
                ciData(trialNum).trialWindowIndex(1) 1]) ...
                max([round(behavData(trialNum).movementTime(1))-...
                ciData(trialNum).trialWindowIndex(1) 1]) + round(behavData(trialNum).movementTime(2))-...
                round(behavData(trialNum).movementTime(1))] ;
            if trialNum==trialId
                % nothing anymore
            end
            %extract raw signal for each trial
            if ~isempty(ciData(trialNum).trialWindowIndex(1):ciData(trialNum).trialWindowIndex(2))
                %             trialWindowIndex=round(behavData(trialNum).movementTime(1))-preAlignWindow:...
                %                 round(behavData(trialNum).movementTime(2))+postAlignWindow;
                
                % convert movement time window to frame time
                timeArray=find(rawTraces{session}.Time>=ciData(trialNum).trialWindowIndex(1),1) :...
                    find(rawTraces{session}.Time<=ciData(trialNum).trialWindowIndex(end),1,'last');
                ciData(trialNum).rawTraceEpochs=getsamples(rawTraces{session},timeArray);
                % get spikes if they come from that trial
                ciData(trialNum).spikes=cellfun(@(eachCell) eachCell(eachCell(:,1)>=ciData(trialNum).trialWindowIndex(1) &...
                    eachCell(:,1)<=ciData(trialNum).trialWindowIndex(end),:), allSpikes, 'UniformOutput', false);
            else
                timeArray=round(behavData(trialNum).movementTime(1)):100:round(behavData(trialNum).movementTime(2));
                ciData(trialNum).rawTraceEpochs=timeseries(NaN(size(timeArray,2),size(rawTraces{1, 1}.Data,2)),timeArray);
                ciData(trialNum).spikes=cell(1,1);
                continue
            end
        end
        trialId=trialNum+1;
    catch
        continue
    end
end

% subtract the neuropil dF/F0 from the somatic dF/F0 to avoid neuropil contamination?
% see Montijn et al. 2016 supp info for example

cd(exportFolder)
save([subject 'Data'],'ciData','rawTraces','calciumEvents','cellIDs','behavData','subject');

% figure; hold on
% foo={ciData.spikes};
% find(cellfun(@(x) ~isempty(x{1}), foo),1)
% plot(ciData(143).rawTraceEpochs.Data(:,1))
% 
% figure; hold on
% plot(MinCIData.F73d1{2, 3}(MinCIData.F73d1{2, 7}(1,1)-200:MinCIData.F73d1{2, 7}(1,2)+200))
% plot(MinCIData.F73d1{2, 3}(MinCIData.F73d1{2, 7}(2,1)-200:MinCIData.F73d1{2, 7}(2,2)+200))
% plot(MinCIData.F73d1{2, 3}(MinCIData.F73d1{2, 7}(3,1)-200:MinCIData.F73d1{2, 7}(3,2)+200))
% figure; hold on
% plot(MinCIData.F73d1{2, 4}(MinCIData.F73d1{2, 7}(1,1)-200:MinCIData.F73d1{2, 7}(1,2)+200))
% plot(MinCIData.F73d1{2, 4}(MinCIData.F73d1{2, 7}(2,1)-200:MinCIData.F73d1{2, 7}(2,2)+200))
% plot(MinCIData.F73d1{2, 4}(MinCIData.F73d1{2, 7}(3,1)-200:MinCIData.F73d1{2, 7}(3,2)+200))



