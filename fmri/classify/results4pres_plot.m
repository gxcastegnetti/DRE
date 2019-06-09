%%
clear
close all

load('results_100permConf_25.mat')

accuracy_between = (acc_BF + acc_FB)/2;
accuracy_within  = (acc_BB + acc_FF)/2;

useSubs = [1:27 29:31];
useSubs = 1:30;
numRois = size(accuracy_between,1);

accuracy_between = accuracy_between(:,useSubs);
accuracy_within = accuracy_within(:,useSubs);

accuracy_within_mean = mean(accuracy_within,2) - 0.5;
accuracy_between_mean = mean(accuracy_between,2) - 0.5;

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
    
    [pCompare(i),~,~] = signrank(accuracy_within(i,:) - accuracy_between(i,:),0,'Tail','right');
end

offset = 0.2;
bar((1:numRois)-offset,accuracy_within_mean,'facecolor',[0.15 0.45 0.75],'barwidth',0.3),hold on
bar((1:numRois)+offset,accuracy_between_mean,'facecolor',[0.55 0.55 0.55],'barwidth',0.3)
errorbar((1:numRois)+offset,accuracy_between_mean,accuracy_between_sem,'k.')
errorbar((1:numRois)-offset,accuracy_within_mean,accuracy_within_sem,'k.')

set(gca,'xtick',1:numRois,'xticklabels',{'vmPFC','OFC','dlPFC'},'fontsize',18)
% title(roiNamesTrue{r},'fontsize',18)
plot(0.5:0.01:50,0.5*ones(length([0.5:0.01:50]),1),'color',[0 0 0],'linestyle','--','linewidth',2)
ylim([-0.01 0.06]),xlim([0.5 3.5]),xtickangle(45)
ylabel('Accuracy above chance')

legend('Within goal','Between goals','location','northwest')

%%
roiNamesTrue = {'vmPFC','OFC','dlPFC'};
myColors = [101,101,101; 101,101,101]/255;
myColors_ss = [201,201,201; 201,201,201]/255;
figure('color',[1 1 1])

hb = bar([accuracy_within_mean,accuracy_between_mean],0.6); hold on
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    
    for jb = 1:numel(roiNames)
        if ib == 1
        scatter(xData(jb)+0.03*randn(numel(subs),1),accuracy_within(jb,:)-0.5,25,'MarkerEdgeColor',myColors_ss(ib,:),...
            'MarkerFaceColor',myColors_ss(ib,:))
        else
           scatter(xData(jb)+0.03*randn(numel(subs),1),accuracy_between(jb,:)-0.5,25,'MarkerEdgeColor',myColors_ss(ib,:),...
            'MarkerFaceColor',myColors_ss(ib,:)) 
        end
    end
    if ib == 1
        errorbar(xData',accuracy_within_mean,accuracy_within_sem,'linestyle','none','color','k','linewidth',2.5,'capsize',0)
    else
        errorbar(xData',accuracy_between_mean,accuracy_between_sem,'linestyle','none','color','k','linewidth',2.5,'capsize',0)
    end
    
    
    hb(ib).FaceColor = myColors(ib,:);
end, clear ib

legend({'Within goal classification','Cross-classification'},'location','northwest'), legend boxoff
set(gca,'fontsize',18,'xtick',1:numel(roiNamesTrue),...
    'xticklabels',roiNamesTrue), ylim([-0.1 0.2]),xtickangle(45)
ylabel('Accuracy (exceedance from chance)')



