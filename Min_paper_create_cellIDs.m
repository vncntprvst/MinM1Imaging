%paste all values
max(foo)
holder=zeros(221,3);foo=[];
% paste Naive / SE1 values
holder(foo(:,2),1:2)=foo;foo=[];
% paste SE1 / SE2 values
holder(foo(:,1),2:3)=foo;foo=[];

holder=holder(logical(sum(holder,2)),:);

cellIndices=table(holder(:,1),holder(:,2),holder(:,3),'VariableNames',{'Early','Expert1','Expert2'});