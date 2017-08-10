function catPct=Min_paper_categorize_neurons(stats)

% responsive neurons, by session
S1TaskRelated=[stats.taskRelated];S1TaskRelated=[S1TaskRelated(1:3:end).indices];
catPct.S1.taskRelated=sum(S1TaskRelated)/numel(S1TaskRelated);

S2TaskRelated=[stats.taskRelated];S2TaskRelated=[S2TaskRelated(2:3:end).indices];
catPct.S2.taskRelated=sum(S2TaskRelated)/numel(S2TaskRelated);

S3TaskRelated=[stats.taskRelated];S3TaskRelated=[S3TaskRelated(3:3:end).indices];
catPct.S3.taskRelated=sum(S3TaskRelated)/numel(S3TaskRelated);


% type of responses
% Session 1
S1categories=[stats.taskRelated]; S1categories=vertcat(S1categories(1:3:end).classification);
earlyCells=cellfun(@(category) contains(category,'early'), S1categories);
preMovementCells=cellfun(@(category) contains(category,'pre'),S1categories); 
movementCells=cellfun(@(category) contains(category,'during'),S1categories);   
lateResponseCells=cellfun(@(category) contains(category,'late'),S1categories);
       
catPct.S1.taskRelated=[sum(earlyCells)/numel(earlyCells), sum(preMovementCells)/numel(preMovementCells),...
    sum(movementCells)/numel(movementCells),sum(lateResponseCells)/numel(lateResponseCells)];
catPct.reachNeurons(1)=(catPct.S1.taskRelated(3))/sum(catPct.S1.taskRelated);
% catPct.reachNeurons(1)=(catPct.S1.taskRelated(2)+catPct.S1.taskRelated(3))/sum(catPct.S1.taskRelated); For Pre-mov + Mov 
catPct.S1.taskRelated=[catPct.S1.taskRelated 1-sum(catPct.S1.taskRelated)];% for complete pie chart

% Session 2
S2categories=[stats.taskRelated];S2categories=vertcat(S2categories(2:3:end).classification);
earlyCells=cellfun(@(category) contains(category,'early'), S2categories);
preMovementCells=cellfun(@(category) contains(category,'pre'), S2categories); 
movementCells=cellfun(@(category) contains(category,'during'), S2categories);   
lateResponseCells=cellfun(@(category) contains(category,'late'), S2categories);
       
catPct.S2.taskRelated=[sum(earlyCells)/numel(earlyCells), sum(preMovementCells)/numel(preMovementCells),...
    sum(movementCells)/numel(movementCells),sum(lateResponseCells)/numel(lateResponseCells)];
catPct.reachNeurons(2)=(catPct.S2.taskRelated(3))/sum(catPct.S2.taskRelated); 
% catPct.reachNeurons(2)=(catPct.S2.taskRelated(2)+catPct.S2.taskRelated(3))/sum(catPct.S2.taskRelated); For Pre-mov + Mov 
catPct.S2.taskRelated=[catPct.S2.taskRelated 1-sum(catPct.S2.taskRelated)];

% Session 3
S3categories=[stats.taskRelated];S3categories=vertcat(S3categories(3:3:end).classification);
earlyCells=cellfun(@(category) contains(category,'early'), S3categories);
preMovementCells=cellfun(@(category) contains(category,'pre'), S3categories); 
movementCells=cellfun(@(category) contains(category,'during'), S3categories);   
lateResponseCells=cellfun(@(category) contains(category,'late'), S3categories);
       
catPct.S3.taskRelated=[sum(earlyCells)/numel(earlyCells), sum(preMovementCells)/numel(preMovementCells),...
    sum(movementCells)/numel(movementCells),sum(lateResponseCells)/numel(lateResponseCells)];
catPct.reachNeurons(3)=(catPct.S3.taskRelated(3))/sum(catPct.S3.taskRelated);
% catPct.reachNeurons(3)=(catPct.S3.taskRelated(2)+catPct.S3.taskRelated(3))/sum(catPct.S3.taskRelated);  For Pre-mov + Mov 
catPct.S3.taskRelated=[catPct.S3.taskRelated 1-sum(catPct.S3.taskRelated)];

% plot proportions
figure;
labels = {'Early','Pre-mov.','Movement','Late response','Non responsive'};
ax1 = subplot(2,3,1);
pie(ax1,catPct.S1.taskRelated,labels)
title(ax1,'Naive session');

ax2 = subplot(2,3,2);
pie(ax2,catPct.S2.taskRelated,labels)
title(ax2,'Expert session 1');

ax3 = subplot(2,3,3);
pie(ax3,catPct.S3.taskRelated,labels)
title(ax3,'Expert session 2');

ax4 = subplot(2,3,5);
plot(catPct.reachNeurons);
set(gca,'XTickLabel',{'NS','ES1','ES2'});
title(ax4,'Proportion of reach neurons (Mov. only) within task related neurons') 

%% categories by outcome-related (success or failure)
% for subject=1:size(stats,2)
%     for session=1:3
%         pctByOutcome(session,subject)=sum([stats(subject).outcomeRelated(session).indices])/...
%             numel([stats(subject).outcomeRelated(session).indices]);
%     end
% end
% catPct.S1.outcomeRelated=mean(pctByOutcome(1,:));
% catPct.S2.outcomeRelated=mean(pctByOutcome(2,:));
% catPct.S3.outcomeRelated=mean(pctByOutcome(3,:));
