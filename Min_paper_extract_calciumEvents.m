%(spike data from Min is already contained in data.ciData.spikes structure) 
% (calcium events foudn in data.calciumEvents )
% [eg: use CaBBI method]
% spikes=CaBBI_method(data(subjectNum).ciData,22.6);
% 
%     for unitNum=1:length(units)  
%     %     spikeTimes=;
%     end

% event detection: 
% three times the median absolute deviation (MAD))
% then non-negative deconvolution method - Huber 2012
% three times the standard deviation - Komiyama 2010 

% Threshold : Min . Deconvolution: OASIS 
% example for neuron #1, session 1
session=1; neuron =1;
exampleData=data(subjectNum).rawTraces{session,1};
trace=exampleData.Data(:,neuron); 
g = 0.92; lambda=2.4; 
[c_oasis_ar1, s_oasis_ar1] = deconvolveCa(trace, 'ar1', g, 'foopsi', 'lambda', lambda); 
%get spikes from Min's threshold detection
cellNum1_spikeTimes=cellfun(@(spikeTimes) spikeTimes(neuron), ...
{data(subjectNum).ciData([data(subjectNum).behavData.session]'==session).spikes},...
'UniformOutput',false); 
cellNum1_spikeTimes=cell2mat([cellNum1_spikeTimes{:}]')/100;

figure;
subplot(2,1,1)
plot(trace); hold on
plot(c_oasis_ar1)
plot(cellNum1_spikeTimes(:,2),ones(size(cellNum1_spikeTimes,1),1),'d') %plot spikes from Min's threshold detection
subplot(2,1,2)
plot(zeros(1,numel(s_oasis_ar1))); hold on
plot(s_oasis_ar1)
plot(cellNum1_spikeTimes(:,2),zeros(size(cellNum1_spikeTimes,1),1),'d') %plot spikes from Min's threshold detection