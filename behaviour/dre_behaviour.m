%% dre_behaviour
% ~~~
% This script:
% 1) Plots subjective rating for each item and condition.
% 2) Plots subjective distribution of ratings in (1).
% 3) Computes and plots logistic regression of choices during day 2 versus
%    the difference between ratings on day 1.
% 4) Plots RDMs based on ratings. These are saved and used as models for
%    the RSA analysis.
% ~~~
% GX Castegnetti 2018

clear
close all
restoredefaultpath

%% folders
fs      = filesep;
dirHere = pwd;
idcs    = strfind(dirHere,fs);
dirNew  = dirHere(1:idcs(end-1)-1); clear mydir idcs RESTOREDEFAULTPATH_EXECUTED
dir.beh = [dirNew,fs,'data',fs,'behaviour'];
dir.psy = [dirNew,fs,'data',fs,'fmri',fs,'psychOut'];

%% subjects and trials
subs = [4:5 8 9 13:17 19:21 23 25:27 29:32 34:37 39];
ntrials = 120;

%% task order (1: fire-boat-fire-boat; 2: boat-fire-boat-fire)
taskOrd = [ones(1,9),2*ones(1,12),1,1,2,1];

%% plots settings
plot_rda_SS = false;
plot_his_SS = 1;
plot_den_SS = false;
hist_fire_color = [0.75 0.5 0];
hist_boat_color = [0.35 0.4 0.8];

%% read tables
objVersion  = 7; % set which column to read according to the object set version
objs        = readtable([dir.beh,fs,'Objects.csv']);
foo         = logical(table2array(objs(:,objVersion)));
objsName    = table2cell(objs(foo,2)); clear foo objs

%% allocate memory
rsmStack_F = nan(length(subs),ntrials,ntrials);
rsmStack_B = nan(length(subs),ntrials,ntrials);
rat_F      = nan(length(subs),ntrials);
rat_B      = nan(length(subs),ntrials);
con_F      = nan(length(subs),ntrials);
con_B      = nan(length(subs),ntrials);
choAll_B   = [];
choAll_F   = [];
ratDiffAll_F = [];
ratDiffAll_B = [];

%% initialise figures
plotVal = figure('color',[1 1 1]);
plotCon = figure('color',[1 1 1]);
plotFam = figure('color',[1 1 1]);
plotCho = figure('color',[1 1 1]);

%% loop over subs
for s = 1:length(subs)
    
    %% day 1 - extract data
    
    % extract data matrices
    data_F = csvread([dirNew,'/data/behaviour/SF',num2str(subs(s),'%03d'),'/SF',num2str(subs(s),'%03d'),'_B1_DRE.csv']); % value/conf. fire
    data_B = csvread([dirNew,'/data/behaviour/SF',num2str(subs(s),'%03d'),'/SF',num2str(subs(s),'%03d'),'_B2_DRE.csv']); % value/conf. boat
    data_p = csvread([dirNew,'/data/behaviour/SF',num2str(subs(s),'%03d'),'/SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']); % familiarity/price
    
    % sort according to object
    [~, idx_F] = sort(data_F(:,2));
    [~, idx_B] = sort(data_B(:,2));
    [~, idx_p] = sort(data_p(:,2));
    
    % extract rating
    rat_F(s,:) = data_F(idx_F,3);
    rat_B(s,:) = data_B(idx_B,3);
    
    % extract confidence
    con_F(s,:) = data_F(idx_F,4);
    con_B(s,:) = data_B(idx_B,4);
    
    % extract familiarity and price
    fam(s,:) = data_p(idx_p,3);
    pri(s,:) = data_p(idx_p,4);
    
    %% plot value
    
    % plot rating for each item
    %     figure('color',[1 1 1])
    %     subplot(2,1,1), bar(rat_F(s,:),'facecolor',hist_fire_color),title(['sub#',num2str(subs(s))])
    %     set(gca,'XTick',1:ntrials,'fontsize',11,'XtickLabel',objsName),xtickangle(45)
    %     subplot(2,1,2), bar(rat_B(s,:),'facecolor',hist_boat_color)
    %     set(gca,'XTick',1:ntrials,'fontsize',11,'XtickLabel',objsName),xtickangle(45)
    
    % histogram value fire
    plotTight = 7;
    figure(plotVal)
    subplot(5,plotTight*5,1+plotTight*(s-1):3+plotTight*(s-1)),histogram(rat_F(s,:),20,'facecolor',hist_fire_color)
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('value')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % histogram value boat
    subplot(5,plotTight*5,4+plotTight*(s-1):6+plotTight*(s-1)),histogram(rat_B(s,:),20,'facecolor',hist_boat_color)
    title(['Sub#',num2str(subs(s)),' - B'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('value')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    
    %% plot confidence
    
    % plot rating for each item
    %     figure('color',[1 1 1])
    %     subplot(2,1,1), bar(con_F(s,:),'facecolor',hist_fire_color),title(['sub#',num2str(subs(s))])
    %     set(gca,'XTick',1:ntrials,'fontsize',11,'XtickLabel',objsName),xtickangle(45)
    %     subplot(2,1,2), bar(con_B(s,:),'facecolor',hist_boat_color)
    %     set(gca,'XTick',1:ntrials,'fontsize',11,'XtickLabel',objsName),xtickangle(45)
    
    % histogram value fire
    figure(plotCon)
    subplot(5,plotTight*5,1+plotTight*(s-1):3+plotTight*(s-1)),histogram(con_F(s,:),20,'facecolor',hist_fire_color)
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('confid.')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % histogram value boat
    subplot(5,plotTight*5,4+plotTight*(s-1):6+plotTight*(s-1)),histogram(con_B(s,:),20,'facecolor',hist_boat_color)
    title(['Sub#',num2str(subs(s)),' - B'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('confid.')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    
    %% plot familiarity
    
    % plot familiarity for each item
    %     figure('color',[1 1 1])
    %     subplot(2,1,1), bar(fam(s,:),'facecolor',hist_fire_color),title(['sub#',num2str(subs(s))])
    %     set(gca,'XTick',1:ntrials,'fontsize',11,'XtickLabel',objsName),xtickangle(45)
    
    % histogram value fire
    figure(plotFam)
    subplot(5,5,s),histogram(fam(s,:),20,'facecolor',[0.25 0.25 0.25])
    title(['Sub#',num2str(subs(s))],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',0:10:50), xlabel('famil.')
    
    %% day 2 - choices
    
    % extract data from the four sessions (1: #trial; 2: TrialStart; 3: TrialType; 4-5: Pic(s)ID; 6: TrialEnd; 7: KeyID; 8: Latency)
    if taskOrd(s) == 1
        typeFirst = 'F';
        typeSecon = 'B';
    else
        typeFirst = 'B';
        typeSecon = 'F';
    end
    
    data_mri_1 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B1',typeFirst,'.csv']);
    data_mri_2 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B2',typeSecon,'.csv']);
    data_mri_3 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B3',typeFirst,'.csv']);
    data_mri_4 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B4',typeSecon,'.csv']);
    
    % concatenate sessions
    if taskOrd(s) == 1
        data_mri_F = [data_mri_1; data_mri_3];
        data_mri_B = [data_mri_2; data_mri_4];
    else
        data_mri_F = [data_mri_2; data_mri_4];
        data_mri_B = [data_mri_1; data_mri_3];
    end
    
    % extract choices
    idxCho_F = data_mri_F(:,3) == 1;
    idxCho_B = data_mri_B(:,3) == 1;
    Cho_F = [data_mri_F(idxCho_F,4:5) data_mri_F(idxCho_F,7)];
    Cho_B = [data_mri_B(idxCho_B,4:5) data_mri_B(idxCho_B,7)];
    
    %% regression rating vs. choice
    
    for i = 1:length(Cho_F)
        
        % extract ID of the objects being displayed - F
        item_id_F_L = Cho_F(i,1);
        item_id_F_R = Cho_F(i,2);
        
        % compute difference between rating - F
        item_rating_F_L = data_F(data_F(:,2) == item_id_F_L,3);
        item_rating_F_R = data_F(data_F(:,2) == item_id_F_R,3);
        ratingDiff_F(i) = item_rating_F_R - item_rating_F_L;
        
        % extract ID of the objects being displayed - B
        item_id_B_L = Cho_B(i,1);
        item_id_B_R = Cho_B(i,2);
        
        % compute difference between rating - B
        item_rating_B_L = data_B(data_B(:,2) == item_id_B_L,3);
        item_rating_B_R = data_B(data_B(:,2) == item_id_B_R,3);
        ratingDiff_B(i) = item_rating_B_R - item_rating_B_L;
    end
    
    % column vectors are nicer
    ratingDiff_F = ratingDiff_F';
    ratingDiff_B = ratingDiff_B';
    
    % logistic regressions
    mdl_F_rht = fitglm(ratingDiff_F,(Cho_F(:,3)+1)/2,'interactions','Distribution','binomial');
    mdl_B_rht = fitglm(ratingDiff_B,(Cho_B(:,3)+1)/2,'interactions','Distribution','binomial');
    mdl_F_opp = fitglm(ratingDiff_F,(Cho_B(:,3)+1)/2,'interactions','Distribution','binomial');
    mdl_B_opp = fitglm(ratingDiff_B,(Cho_F(:,3)+1)/2,'interactions','Distribution','binomial');
    
    
    %% plot logistic curves of choice vs dV
    
    % draw sigmoids with the fitted parameters
    xspan = -50:0.1:50;
    sigm_F_rht = sigmf(xspan,[mdl_F_rht.Coefficients.Estimate(2) mdl_F_rht.Coefficients.Estimate(1)]);
    sigm_B_rht = sigmf(xspan,[mdl_B_rht.Coefficients.Estimate(2) mdl_B_rht.Coefficients.Estimate(1)]);
    sigm_F_opp = sigmf(xspan,[mdl_F_opp.Coefficients.Estimate(2) mdl_F_opp.Coefficients.Estimate(1)]);
    sigm_B_opp = sigmf(xspan,[mdl_B_opp.Coefficients.Estimate(2) mdl_B_opp.Coefficients.Estimate(1)]);
    
    % choice during fire
    figure(plotCho)
    subplot(5,plotTight*5,1+plotTight*(s-1):3+plotTight*(s-1))
    plot(xspan,sigm_F_rht,'linewidth',5,'color',hist_fire_color), hold on % based on value assigned during fire
    plot(xspan,sigm_F_opp,'linewidth',5,'color',hist_boat_color) % based on value assigned during boat
    plot(ratingDiff_F,(Cho_F(:,3)+1)/2,'linestyle','none','marker','.','markersize',15,'color','k')
    set(gca,'fontsize',12,'ytick',[],'xtick',[-50 0 50]),ylabel('P(Right)'),xlabel('Value diff.')
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14)
%     legend('F-based','B-based')
        
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.25*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % choice during fire
    subplot(5,plotTight*5,4+plotTight*(s-1):6+plotTight*(s-1))
    plot(xspan,sigm_B_rht,'linewidth',5,'color',hist_boat_color),hold on % based on value assigned during boat
    plot(xspan,sigm_B_opp,'linewidth',5,'color',hist_fire_color) % based on value assigned during fire
    plot(ratingDiff_B,(Cho_B(:,3)+1)/2,'linestyle','none','marker','.','markersize',25,'color','k')
    set(gca,'fontsize',12,'ytick',[],'xtick',[-50 0 50]),xlabel('Value diff.')
    title(['Sub#',num2str(subs(s)),' - B'],'fontsize',14) 
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.25*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    
    %% compute RSA (sort of, it's just the abs difference between ratings)
    RSM_fire_ss = nan(ntrials);
    RSM_door_ss = nan(ntrials);
    for i = 1:ntrials
        RSM_fire_ss(:,i) = abs(rat_F(s,i)/50 - rat_F(s,:)/50);
        RSM_door_ss(:,i) = abs(rat_B(s,i)/50 - rat_B(s,:)/50);
    end
    rsmStack_F(s,:,:) = RSM_fire_ss;
    rsmStack_B(s,:,:) = RSM_door_ss;
    
    %% fill structure for RSA toolbox
    %     RDM_fire_struct{s}.RDM = RSM_fire_ss;
    %     RDM_door_struct{s}.RDM = RSM_door_ss;
    %     RDM_fire_struct{s}.name = ['sub#',num2str(subs(s)),'_fire'];
    %     RDM_door_struct{s}.name = ['sub#',num2str(subs(s)),'_door'];
    %     RDM_fire_struct{s}.color = [0 0 1];
    %     RDM_door_struct{s}.color = [0 0 1];
    %
    %     if pl_rda_SS
    %         figure(hRSA)
    %         subplot(2,length(subs),s), imagesc(RSM_fire_ss)%,caxis([0 4])
    %         subplot(2,length(subs),s+length(subs)), imagesc(RSM_door_ss)%,caxis([0 4])
    %     end
    %
    %     clear ratingsRDA ratings_z confids_z
    %
    %     % compute correlation between ratings and confidences (include between correlation and confidence)
    %     fooRat = ~isnan(rat_F) & ~isnan(rat_B);
    %     fooCon = ~isnan(con_F) & ~isnan(con_B);
    %     ratings_fire_ = rat_F(fooRat);
    %     ratings_door_ = rat_B(fooRat);
    %     confids_fire_ = con_F(fooCon);
    %     confids_door_ = con_B(fooCon);
    %     [r_rat,p_rat] = corrcoef(ratings_fire_,ratings_door_);
    %     [r_con,p_con] = corrcoef(confids_fire_,confids_door_);
    %     R_rating(s) = r_rat(2,1);
    %     p_rating(s) = p_rat(2,1);
    %     R_confid(s) = r_con(2,1);
    %     p_confid(s) = p_con(2,1);
    clear r_rat p_rat r_con p_con ratings_fire_ ratings_door_ confids_fire_ confids_door_
    
end
keyboard
%% table with single subject results
TableCorr = table(subs',R_rating',p_rating',R_confid',p_confid');
TableCorr.Properties.VariableNames = {'Subject' 'R_rating' 'p_rating' 'R_confid','p_confid'};
clear R_rating p_rating R_confid p_confid rda_SS hist_SS i blk rda_SS hist_SS zzScore

%% table with single subject results
TableCorr = table(subs',rViv_F',pViv_F',rViv_B',pViv_B');
TableCorr.Properties.VariableNames = {'Subject' 'rVivFire' 'pVivFire' 'rVivDoor','pVivDoor'};
clear R_rating p_rating R_confid p_confid rda_SS hist_SS i blk rda_SS hist_SS zzScore

%% compute pairwise correlation
ratRDA_all = cat(1,rsmStack_F,rsmStack_B);
j = 1;
for s1 = 1:2*length(subs)
    if j < length(subs) + 1
        sessNames{j} = ['sub#',num2str(subs(j))];
        sessNames{j+length(subs)} = ['sub#',num2str(subs(j))];
        j = j + 1;
    end
    for s2 = 1:2*length(subs)
        rda1 = squeeze(ratRDA_all(s1,:,:));
        rda2 = squeeze(ratRDA_all(s2,:,:));
        fooNaN = isnan(rda1) | isnan(rda2);
        [r,p] = corrcoef(rda1(~fooNaN),rda2(~fooNaN));
        r_xxx(s1,s2) = r(2,1); %#ok<*SAGROW>
        p_xxx(s1,s2) = p(2,1);
    end
end
figure('color',[1 1 1])
imagesc(r_xxx),caxis([0 1])
set(gca,'XTick',1:2*length(subs),'YTick',1:2*length(subs),'fontsize',12,...
    'XtickLabel',sessNames,'YtickLabel',sessNames), xtickangle(45)
clear s1 s2 j sessNames rda1 rda2 fooNaN r p

%% plot RDA
ratingsFire_mean = (rat_F);
[~,idx] = sort(rat_F,'descend');
ratingsDoor_mean = nanmean(rat_B);
ratRDA_fire_mean = squeeze(nanmean(rsmStack_F,1));
ratRDA_door_mean = squeeze(nanmean(rsmStack_B,1));
figure('color',[1 1 1])
imagesc(RSM_fire_ss),set(gca,'XTick',1:100,'YTick',1:100,'fontsize',10,...
    'XtickLabel',objsName,'YtickLabel',objsName)
xtickangle(45),caxis([0 1])
figure('color',[1 1 1])
imagesc(RSM_door_ss),set(gca,'XTick',1:100,'YTick',1:100,'fontsize',10,...
    'XtickLabel',objsName,'YtickLabel',objsName)
xtickangle(45),caxis([0 1])

%% apply real RDA (from toolbox)
userOptions.analysisName = 'DRE_test'; % string prepended to the saved files
userOptions.rootPath = '/Users/gcastegnetti/Dropbox/DRE/analysis'; % string describing the root path where files will be
userOptions.rankTransform = false;
userOptions.conditionLabels = objsName;
if plot_den_SS
    for s = 1:length(subs)
        dendrogramConditions(RDM_fire_struct{s}, userOptions)
        dendrogramConditions(RDM_door_struct{s}, userOptions)
    end
end

% plot dendogram of averaged RDMs
meanRatRDM_fire.RDM = ratRDA_fire_mean;
meanRatRDM_fire.name = 'Fire';
meanRatRDM_fire.color = [0 0 1];

meanRatRDM_door.RDM = ratRDA_door_mean;
meanRatRDM_door.name = 'Boat';
meanRatRDM_door.color = [0 0 1];

dendrogramConditions(meanRatRDM_fire, userOptions)
dendrogramConditions(meanRatRDM_door, userOptions)
