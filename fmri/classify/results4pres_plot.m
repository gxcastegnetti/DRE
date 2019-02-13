%%
clear
close all

load('results_1000perm_33.mat')

accuracy_between = (acc_BF + acc_FB)/2;
accuracy_within  = (acc_BB + acc_FF)/2;

useSubs = [1:27 29:31];
numRois = size(accuracy_between,1);

accuracy_between = accuracy_between(:,useSubs);
accuracy_within = accuracy_within(:,useSubs);

accuracy_within_mean = mean(accuracy_within,2);
accuracy_between_mean = mean(accuracy_between,2);

accuracy_within_sem = std(accuracy_within,0,2)/sqrt(numel(useSubs));
accuracy_between_sem = std(accuracy_between,0,2)/sqrt(numel(useSubs));

figure('color',[1 1 1])
for i = 1:numRois
    
    % within context
%     [~,pWithin(i),~,~] = ttest(accuracy_within(i,:)-0.5,0,'Tail','right');
    [pWithin(i),~,~] = signrank(accuracy_within(i,:)-0.5,0,'Tail','right');
    
    % between context
%     [~,pBetween(i),~,~] = ttest(accuracy_between(i,:)-0.5,0,'Tail','right');
    [pBetween(i),~,~] = signrank(accuracy_between(i,:)-0.5,0,'Tail','right');
end

offset = 0.2;
bar((1:numRois)-offset,accuracy_within_mean,'facecolor',[0.15 0.45 0.75],'barwidth',0.3),hold on
bar((1:numRois)+offset,accuracy_between_mean,'facecolor',[0.55 0.55 0.55],'barwidth',0.3)
errorbar((1:numRois)+offset,accuracy_between_mean,accuracy_between_sem,'k.')
errorbar((1:numRois)-offset,accuracy_within_mean,accuracy_within_sem,'k.')

set(gca,'xtick',1:numRois,'xticklabels',{'vmPFC','OFC','dlPFC'},'fontsize',18)
% title(roiNamesTrue{r},'fontsize',18)
plot(0.5:0.01:50,0.5*ones(length([0.5:0.01:50]),1),'color',[0 0 0],'linestyle','--','linewidth',2)
ylim([0.46 0.58]),xlim([0.5 3.5]),xtickangle(45)
ylabel('Classification accuracy')

legend('Within goal','Between goals')

%% symmetry index
% simmIdx = -accuracy_between + accuracy_within;
% simmIdx_mean = mean(simmIdx,2);
% simmIdx_sem = std(simmIdx,0,2)/sqrt(numel(useSubs));
% 
% for i = 1:numRois
%     [pDiff(i),~,~] = signrank(simmIdx(i,:),0,'Tail','right');
% end

%% symmetry index
simmIdx = accuracy_within;
simmIdx_mean = mean(simmIdx,2);
simmIdx_sem = std(simmIdx,0,2)/sqrt(numel(useSubs));

for i = 1:numRois
    [pDiff(i),~,~] = signrank(simmIdx(i,:),0,'Tail','right');
end


