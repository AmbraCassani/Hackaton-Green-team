function [] = plot_Boxplot(struct1,struct2,field)
% Boxplots

% make a boxplot for x:
figure
boxplot(struct1.(field))
hold on
plot(ones(length(struct1.(field)),1), struct1.(field),'O')
% make a boxplot for y, and specify that it belongs at 2 along the x-axis:
hold on
boxplot(struct2.(field),'Positions',2)
hold on
plot(ones(length(struct2.(field)),1)*2, struct2.(field),'O')
% update axes x-ticks and labels:
set(gca(),'XTick',[1 2],'XTickLabels',{'Control','Fallers'})
end