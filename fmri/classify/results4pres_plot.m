load('results4pres2.mat')

useRois = [2 4 5 9 10];
useRois = 1:8;

accuracy_between = accuracy_between(useRois,:);
accuracy_within = accuracy_within(useRois,:);

accuracy_within_mean = mean(accuracy_within,2);
accuracy_between_mean = mean(accuracy_between,2);

accuracy_within_sem = std(accuracy_within,0,2)/sqrt(30);
accuracy_between_sem = std(accuracy_between,0,2)/sqrt(30);

figure('color',[1 1 1])
for i = 1:numel(useRois)
    
    % within context
    [~,pWithin(i),~,~] = ttest(accuracy_within(i,:)-0.5,0,'Tail','right');
    
    % between context
    [~,pBetween(i),~,~] = ttest(accuracy_between(i,:)-0.5,0,'Tail','right');
    
end

bar(1:8,accuracy_within_mean,'facecolor',[0.15 0.45 0.75],'barwidth',0.3),hold on
errorbar(1:8,accuracy_within_mean,accuracy_within_sem,'k.')
    
bar(1.36:8.36,accuracy_between_mean,'facecolor',[0.55 0.55 0.55],'barwidth',0.3)
errorbar(1.36:8.36,accuracy_between_mean,accuracy_between_sem,'k.')

set(gca,'xtick',1.18:8.18,'xticklabels',{'Lingual','lHPC','rHPC','ACC','rIns','vmPFC','pOFC','aOFC'},'fontsize',18)
% title(roiNamesTrue{r},'fontsize',18)
plot(0.5:0.01:50,0.5*ones(length([0.5:0.01:50]),1),'color',[0 0 0],'linestyle','--','linewidth',2)
ylim([0.46 0.58]),xlim([0.5 8.8]),xtickangle(45)
ylabel('Classification accuracy')
