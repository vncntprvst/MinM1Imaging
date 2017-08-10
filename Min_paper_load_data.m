function [data,averages,stats]=Min_paper_load_data(preAlignWindow,postAlignWindow,samplingRate,subjectList)

cd ('D:\Sync\Box Sync\Home Folder vp35\Sync\other research\Min\GCamp6_movement')

%% load data
addpath(pwd);
dataFolder = [pwd filesep 'Data' filesep];
addpath(dataFolder);

if exist([dataFolder 'all_data.mat'],'file')
    load([dataFolder 'all_data.mat']);
else
    
    data=struct('ciData',[],'rawTraces',{},'calciumEvents',[],...
        'cellIDs',[],'behavData',[],'subject',[]);
    numSubject=size(subjectList,1);
    for subjectNum=1:numSubject
        % if file doesn't exist, import and create file first
        if ~exist([dataFolder subjectList{subjectNum} 'Data.mat'],'file')
            Min_paper_getCIdata(dataFolder,[preAlignWindow postAlignWindow],samplingRate,subjectList{subjectNum});
        end
        data(subjectNum)=...
            load([dataFolder subjectList{subjectNum} 'Data.mat'],...
            'ciData','rawTraces','calciumEvents','cellIDs',...
            'behavData','subject');
    end
    [averages,stats]=deal([]);
end
