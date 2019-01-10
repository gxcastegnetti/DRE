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
subs = [4:5 7 8 9 13:17 19:21 23 25:26 29:32 34:35 37 39:41 43 47:50];
% subs = [4:5 7 8 9 10 13:17 19 20 21 23 25:26 27 29:32 34:35 37 39:41 43 44 45 47:49 50];
ntrials = 120;

%% task order (1: fire-boat-fire-boat; 2: boat-fire-boat-fire)
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% plots settings
plot_rda_SS = false;
plot_his_SS = 1;
plot_den_SS = false;
hist_fire_color = [0.8706 0.4196 0.1569];
hist_boat_color = [0.2627 0.5255 0.5882];

%% read tables
objVersion  = 7; % set which column to read according to the object set version
objs        = readtable([dir.beh,fs,'Objects.csv']);
foo         = logical(table2array(objs(:,objVersion)));
objsName    = table2cell(objs(foo,2)); clear foo objs

%% allocate memory
rsmStack_F = nan(length(subs),ntrials,ntrials);
rsmStack_B = nan(length(subs),ntrials,ntrials);
val_F      = nan(length(subs),ntrials);
val_B      = nan(length(subs),ntrials);
con_F      = nan(length(subs),ntrials);
con_B      = nan(length(subs),ntrials);
fam        = nan(length(subs),ntrials);
pri        = nan(length(subs),ntrials);
allSubs    = [];
allCond    = [];
alldV      = [];
alldV_opp  = [];
alldF      = [];
allCho     = [];

%% initialise figures
plotVal = figure('color',[1 1 1]);
plotCon = figure('color',[1 1 1]);
plotFam = figure('color',[1 1 1]);
plotCho = figure('color',[1 1 1]);
plotCrr = figure('color',[1 1 1]);

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
    val_F(s,:) = data_F(idx_F,3);
    val_B(s,:) = data_B(idx_B,3);
    
    % extract confidence
    con_F(s,:) = data_F(idx_F,4);
    con_B(s,:) = data_B(idx_B,4);
    
    % extract familiarity and price
    fam(s,:) = data_p(idx_p,3);
    pri(s,:) = data_p(idx_p,4);
    
    % create matrix useful later for computing correlations
    all{s} = [val_F(s,:)' val_B(s,:)' con_F(s,:)' con_B(s,:)' fam(s,:)'];
    
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
    subplot(6,plotTight*6,1+plotTight*(s-1):3+plotTight*(s-1)),histogram(val_F(s,:),20,'facecolor',hist_fire_color)
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('value')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % histogram value boat
    subplot(6,plotTight*6,4+plotTight*(s-1):6+plotTight*(s-1)),histogram(val_B(s,:),20,'facecolor',hist_boat_color)
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
    subplot(6,plotTight*6,1+plotTight*(s-1):3+plotTight*(s-1)),histogram(con_F(s,:),20,'facecolor',hist_fire_color)
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',[1 50]), xlabel('confid.')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.15*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % histogram value boat
    subplot(6,plotTight*6,4+plotTight*(s-1):6+plotTight*(s-1)),histogram(con_B(s,:),20,'facecolor',hist_boat_color)
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
    subplot(6,6,s),histogram(fam(s,:),20,'facecolor',[0.25 0.25 0.25])
    title(['Sub#',num2str(subs(s))],'fontsize',14), set(gca,'fontsize',12,'ytick',[],'xtick',0:10:50), xlabel('famil.')
    
    
    %% day 2 - extract choices
    
    % extract data from the four sessions (1: #trial; 2: TrialStart; 3: TrialType; 4-5: Pic(s)ID; 6: TrialEnd; 7: KeyID; 8: Latency)
    if taskOrd(s) == 1
        goal1 = 'F';
        goal2 = 'B';
    else
        goal1 = 'B';
        goal2 = 'F';
    end
    
    % extract psychopy data from the four scanning sessions
    data_mri_1 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B1',goal1,'.csv']);
    data_mri_2 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B2',goal2,'.csv']);
    data_mri_3 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B3',goal1,'.csv']);
    data_mri_4 = csvread([dir.psy,fs,'SF',num2str(subs(s),'%03d'),'/DRE_mri_S',num2str(subs(s),'%03d'),'_B4',goal2,'.csv']);
    
    % concatenate sessions
    if taskOrd(s) == 1
        data_mri_F = [data_mri_1; data_mri_3];
        data_mri_B = [data_mri_2; data_mri_4];
    else
        data_mri_F = [data_mri_2; data_mri_4];
        data_mri_B = [data_mri_1; data_mri_3];
    end
    
    % extract indices of choice trials
    idxCho_F = data_mri_F(:,3) == 1;
    idxCho_B = data_mri_B(:,3) == 1;
    
    % vector with 3 columns: L obj, R obj, choice
    choice_F = [data_mri_F(idxCho_F,4:5) data_mri_F(idxCho_F,7:8)];
    choice_B = [data_mri_B(idxCho_B,4:5) data_mri_B(idxCho_B,7:8)];
    
    clear data_mri_1 data_mri_2 data_mri_3 data_mri_4 idxCho_F idxCho_B goal 1 goal2
    
    %% val, con, fam vs RT
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find val, con, fam of the objects %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(choice_F)
        
        % --- find when they were presented on day 1 ---
        
        % fire
        idxCho_F_L(i) = find(data_F(:,2) == choice_F(i,1));
        idxCho_F_R(i) = find(data_F(:,2) == choice_F(i,2));
        
        % boat
        idxCho_B_L(i) = find(data_B(:,2) == choice_B(i,1));
        idxCho_B_R(i) = find(data_B(:,2) == choice_B(i,2));
        
        % --- find val, fam, con assigned on day 1 ---
        
        % fire
        if choice_F(i,3) == -1
            valCho_F(i) = data_F(idxCho_F_L(i),3);
            valUnc_F(i) = data_F(idxCho_F_R(i),3);
            conCho_F(i) = data_F(idxCho_F_L(i),4);
            conUnc_F(i) = data_F(idxCho_F_R(i),4);
            famCho_F(i) = data_p(idxCho_F_L(i),3);
            famUnc_F(i) = data_p(idxCho_F_R(i),3);
        elseif choice_F(i,3) == 1
            valCho_F(i) = data_F(idxCho_F_R(i),3);
            valUnc_F(i) = data_F(idxCho_F_L(i),3);
            conCho_F(i) = data_F(idxCho_F_R(i),4);
            conUnc_F(i) = data_F(idxCho_F_L(i),4);
            famCho_F(i) = data_p(idxCho_F_R(i),3);
            famUnc_F(i) = data_p(idxCho_F_L(i),3);
        end
        
        % boat
        if choice_B(i,3) == -1
            valCho_B(i) = data_B(idxCho_B_L(i),3);
            valUnc_B(i) = data_B(idxCho_B_R(i),3);
            conCho_B(i) = data_B(idxCho_B_L(i),4);
            conUnc_B(i) = data_B(idxCho_B_R(i),4);
            famCho_B(i) = data_p(idxCho_B_L(i),3);
            famUnc_B(i) = data_p(idxCho_B_R(i),3);
        elseif choice_B(i,3) == 1
            valCho_B(i) = data_B(idxCho_B_R(i),3);
            valUnc_B(i) = data_B(idxCho_B_L(i),3);
            conCho_B(i) = data_B(idxCho_B_R(i),4);
            conUnc_B(i) = data_B(idxCho_B_L(i),4);
            famCho_B(i) = data_p(idxCho_B_R(i),3);
            famUnc_B(i) = data_p(idxCho_B_L(i),3);
        end
        
    end
    
    % put conditions together
    valCho = [valCho_F'; valCho_B'];
    valUnc = [valUnc_F'; valUnc_B'];
    conCho = [conCho_F'; conCho_B'];
    conUnc = [conUnc_F'; conUnc_B'];
    famCho = [famCho_F'; famCho_B'];
    famUnc = [famUnc_F'; famUnc_B'];
    
    % take absolute difference
    valDif = abs(valCho - valUnc);
    conDif = abs(conCho - conUnc);
    famDif = abs(famCho - famUnc);
    
    % take chosen - unchosen
    valChMUnc = valCho - valUnc;
    conChMUnc = conCho - conUnc;
    famChMUnc = famCho - famUnc;
    
    %%%%%%%%%%%%
    % take RTs %
    %%%%%%%%%%%%
    
    rt_F = choice_F(:,4);
    rt_B = choice_B(:,4);
    rt = [rt_F;rt_B];
    
    % take means
    rt_F_mean(s) = nanmean(rt_F);
    rt_B_mean(s) = nanmean(rt_B);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % correlations scores/RTs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r_valDif = corr(rt,valDif,'rows','complete');
    r_conDif = corr(rt,conDif,'rows','complete');
    r_famDif = corr(rt,famDif,'rows','complete');
    r_valChMUnc = corr(rt,valChMUnc,'rows','complete');
    r_conChMUnc = corr(rt,conChMUnc,'rows','complete');
    r_famChMUnc = corr(rt,famChMUnc,'rows','complete');
    
    fooPlot_ChMUnc(s,:) = [r_valChMUnc r_conChMUnc r_famChMUnc];
    fooPlot_Chosen(s,:) = [r_valChMUnc r_conChMUnc r_famChMUnc];
    
    % plot
    figure(plotCrr)
    subplot(6,6,s)
    bar(fooPlot_Chosen(s,:))
    set(gca,'fontsize',12,'xticklabel',{'v','c','f'}),ylabel('r')
    title(['Sub#',num2str(subs(s))],'fontsize',14)
    
    
    %% logistic regression of choice on day 2 from rating on day 1
    
    % loop over (48) choice trials for each goal
    for i = 1:length(choice_F)
        
        %%%%%%%%
        % fire %
        %%%%%%%%
        
        % extract ID of the two objects being displayed
        item_id_F_L = choice_F(i,1);
        item_id_F_R = choice_F(i,2);
        
        % dV
        item_rating_F_L = data_F(data_F(:,2) == item_id_F_L,3);
        item_rating_F_R = data_F(data_F(:,2) == item_id_F_R,3);
        dV_F(i) = item_rating_F_R - item_rating_F_L; %#ok<*SAGROW>
        
        % dF
        item_famili_F_L = data_p(data_p(:,2) == item_id_F_L,3);
        item_famili_F_R = data_p(data_p(:,2) == item_id_F_R,3);
        dF_F(i) = item_famili_F_R - item_famili_F_L;
        
        %%%%%%%%
        % boat %
        %%%%%%%%
        
        % extract ID of the two objects being displayed
        item_id_B_L = choice_B(i,1);
        item_id_B_R = choice_B(i,2);
        
        % dV
        item_rating_B_L = data_B(data_B(:,2) == item_id_B_L,3);
        item_rating_B_R = data_B(data_B(:,2) == item_id_B_R,3);
        dV_B(i) = item_rating_B_R - item_rating_B_L;
        
        % dF
        item_famili_B_L = data_p(data_p(:,2) == item_id_B_L,3);
        item_famili_B_R = data_p(data_p(:,2) == item_id_B_R,3);
        dF_B(i) = item_famili_B_R - item_famili_B_L;
        
    end, clear item_id_F_L item_id_F_R item_id_B_L item_id_B_R item_rating_F_L item_rating_F_R item_rating_B_L item_rating_B_R
    clear item_famili_B_L item_famili_B_R item_famili_F_L item_famili_F_R
    
    % column vectors are nicer
    dV_F = dV_F(:);
    dV_B = dV_B(:);
    dF_F = dF_F(:);
    dF_B = dF_B(:);
    
    % logistic regressions
    mdl_FonF = fitglm(dV_F,(choice_F(:,3)+1)/2,'interactions','Distribution','binomial'); % choice during F vs rating during F
    mdl_BonB = fitglm(dV_B,(choice_B(:,3)+1)/2,'interactions','Distribution','binomial'); % choice during B vs rating during B
    mdl_BonF = fitglm(dV_F,(choice_B(:,3)+1)/2,'interactions','Distribution','binomial'); % choice during B vs rating during F
    mdl_FonB = fitglm(dV_B,(choice_F(:,3)+1)/2,'interactions','Distribution','binomial'); % choice during F vs rating during B
    
    % extract p values (just better interpretation at this stage)
    pVal(s).('FonF') = mdl_FonF.Coefficients.pValue(2);
    pVal(s).('FonB') = mdl_FonB.Coefficients.pValue(2);
    pVal(s).('BonB') = mdl_BonB.Coefficients.pValue(2);
    pVal(s).('BonF') = mdl_BonF.Coefficients.pValue(2);
    
    % extract slopes
    slopes(s).('FonF') = mdl_FonF.Coefficients.Estimate(2);
    slopes(s).('FonB') = mdl_FonB.Coefficients.Estimate(2);
    slopes(s).('BonB') = mdl_BonB.Coefficients.Estimate(2);
    slopes(s).('BonF') = mdl_BonF.Coefficients.Estimate(2);
    
    % extract offsets
    offset(s).('FonF') = mdl_FonF.Coefficients.Estimate(1);
    offset(s).('FonB') = mdl_FonB.Coefficients.Estimate(1);
    offset(s).('BonB') = mdl_BonB.Coefficients.Estimate(1);
    offset(s).('BonF') = mdl_BonF.Coefficients.Estimate(1);
    
    %% create vectors for table
    allSubs = [allSubs; s*ones(2*length(choice_F),1)];
    allCond = [allCond; ones(length(choice_F),1); 2*ones(length(choice_F),1)];
    alldV   = [alldV; dV_F/50; dV_B/50];
    alldV_opp   = [alldV_opp; dV_B/50; dV_F/50];
    alldF   = [alldF; dF_F/50; dF_B/50]; %#ok<*AGROW>
    allCho  = [allCho; choice_F(:,3); choice_B(:,3)];
    
    %% plot logistic curves of choice vs dV
    
    % draw sigmoids with the fitted parameters
    xspan = -50:0.1:50;
    sigm_FonF = sigmf(xspan,[mdl_FonF.Coefficients.Estimate(2) mdl_FonF.Coefficients.Estimate(1)]);
    sigm_BonB = sigmf(xspan,[mdl_BonB.Coefficients.Estimate(2) mdl_BonB.Coefficients.Estimate(1)]);
    sigm_BonF = sigmf(xspan,[mdl_BonF.Coefficients.Estimate(2) mdl_BonF.Coefficients.Estimate(1)]);
    sigm_FonB = sigmf(xspan,[mdl_FonB.Coefficients.Estimate(2) mdl_FonB.Coefficients.Estimate(1)]);
    clear mdl_FonF mdl_BonB mdl_FonB mdl_BonF
    
    % choice during fire
    figure(plotCho)
    subplot(6,plotTight*6,1+plotTight*(s-1):3+plotTight*(s-1))
    plot(xspan,sigm_FonF,'linewidth',5,'color',hist_fire_color), hold on % based on value assigned during fire
    plot(xspan,sigm_FonB,'linewidth',5,'color',hist_boat_color) % based on value assigned during boat
    plot(dV_F,(choice_F(:,3)+1)/2,'linestyle','none','marker','.','markersize',15,'color','k')
    set(gca,'fontsize',12,'ytick',[],'xtick',[-50 0 50]),ylabel('P(Right)'),xlabel('Value diff.')
    title(['Sub#',num2str(subs(s)),' - F'],'fontsize',14)
    %     legend('F-based','B-based')
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.25*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    % choice during boat
    subplot(6,plotTight*6,4+plotTight*(s-1):6+plotTight*(s-1))
    plot(xspan,sigm_BonB,'linewidth',5,'color',hist_boat_color),hold on % based on value assigned during boat
    plot(xspan,sigm_BonF,'linewidth',5,'color',hist_fire_color) % based on value assigned during fire
    plot(dV_B,(choice_B(:,3)+1)/2,'linestyle','none','marker','.','markersize',25,'color','k')
    set(gca,'fontsize',12,'ytick',[],'xtick',[-50 0 50]),xlabel('Value diff.')
    title(['Sub#',num2str(subs(s)),' - B'],'fontsize',14)
    
    % make label closer to axis
    xh = get(gca,'xlabel'); % handle to the label object
    p = get(xh,'position'); % get the current position property
    p(2) = 0.25*p(2) ;      % reduce distance,
    set(xh,'position',p)    % set the new position
    
    
    %% compute RSA (sort of, it's just the abs difference between ratings)
    RSM_ss_F = nan(ntrials);
    RSM_ss_B = nan(ntrials);
    for i = 1:ntrials
        RSM_ss_F(:,i) = abs(val_F(s,i)/50 - val_F(s,:)/50);
        RSM_ss_B(:,i) = abs(val_B(s,i)/50 - val_B(s,:)/50);
    end
    rsmStack_F(s,:,:) = RSM_ss_F;
    rsmStack_B(s,:,:) = RSM_ss_B;
    
    figure('color',[1 1 1])
    subplot(1,2,1),imagesc(RSM_ss_F),set(gca,'fontsize',12,'ytick',[],'xtick',[])
    subplot(1,2,2),imagesc(RSM_ss_B),set(gca,'fontsize',12,'ytick',[],'xtick',[])
    colormap bone
    
    %% fill structure for RSA toolbox
    %         RDM_fire_struct{s}.RDM = RSM_fire_ss;
    %         RDM_door_struct{s}.RDM = RSM_door_ss;
    %         RDM_fire_struct{s}.name = ['sub#',num2str(subs(s)),'_F'];
    %         RDM_door_struct{s}.name = ['sub#',num2str(subs(s)),'_B'];
    %         RDM_fire_struct{s}.color = [0 0 1];
    %         RDM_door_struct{s}.color = [0 0 1];
    %
    %         if pl_rda_SS
    %             figure(hRSA)
    %             subplot(2,length(subs),s), imagesc(RSM_fire_ss)%,caxis([0 4])
    %             subplot(2,length(subs),s+length(subs)), imagesc(RSM_door_ss)%,caxis([0 4])
    %         end
    %
    %     clear ratingsRDA ratings_z confids_z
    %
    
    clear r_rat p_rat r_con p_con ratings_fire_ ratings_door_ confids_fire_ confids_door_
    
end

%% correlation between value, confidence, familiarity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% value and confidence between subjects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute correlation matrices
c_val_F = corr(val_F','rows','complete');
c_val_B = corr(val_B','rows','complete');
c_con_F = corr(con_F','rows','complete');
c_con_B = corr(con_B','rows','complete');

% plot
figure('color',[1 1 1])
subplot(2,2,1),imagesc(c_val_F), title('value - F','fontsize',14), ylabel('# sub'), xlabel('# sub'), caxis([-1 1]), set(gca,'fontsize',16), colorbar
subplot(2,2,2),imagesc(c_con_F), title('confid. - F','fontsize',14), ylabel('# sub'), xlabel('# sub'), caxis([-1 1]), set(gca,'fontsize',16), colorbar
subplot(2,2,3),imagesc(c_val_B), title('value - B','fontsize',14), ylabel('# sub'), xlabel('# sub'), caxis([-1 1]), set(gca,'fontsize',16), colorbar
subplot(2,2,4),imagesc(c_con_B), title('confid. - B','fontsize',14), ylabel('# sub'), xlabel('# sub'), caxis([-1 1]), set(gca,'fontsize',16), colorbar

% compute familiarity between subjects
c_fam = corr(fam');

% plot
figure('color',[1 1 1])
imagesc(c_fam), title('familiarity','fontsize',14), ylabel('# sub'), xlabel('# sub'), caxis([-1 1]), set(gca,'fontsize',16), colorbar

clear c_val_F c_val_B c_con_F c_con_B c_fam

%%%%%%%%%%%%%%%%%%%
% within subjects %
%%%%%%%%%%%%%%%%%%%

figure('color',[1 1 1])
for s = 1:length(subs)
    
    % compute correlation between val_f, val_B, con_F, con_B, fam, for each subject
    cMat(:,:,s) = corr(all{s},'rows','complete');
    
    % plot
    subplot(6,6,s)
    imagesc(squeeze(cMat(:,:,s))), title(['sub#',num2str(subs(s))]), caxis([-1 1])
    set(gca,'XTick',1:5,'YTick',1:5,'fontsize',11,'XtickLabel',{'vF','vB','cF','cB','fm'},'YtickLabel',{'vF','vB','cF','cB','fm'})
end, clear s

% plot average
figure('color',[1 1 1])
imagesc(squeeze(mean(cMat,3))), title('average','fontsize',16), caxis([-1 1])
set(gca,'XTick',1:5,'YTick',1:5,'fontsize',16,'XtickLabel',{'vF','vB','cF','cB','fm'},'YtickLabel',{'vF','vB','cF','cB','fm'})

clear cMat all

%% plot average correlation between scores and RT
figure
bar(mean(fooPlot_ChMUnc,1))
set(gca,'fontsize',12,'xticklabel',{'v','c','f'}),ylabel('r')
clear fooPlot_ChMUnc

%% plot average decision sigmoid
slope_avg.FonF = mean([slopes.FonF]);
slope_avg.FonB = mean([slopes.FonB]);
slope_avg.BonF = mean([slopes.BonF]);
slope_avg.BonB = mean([slopes.BonB]);

offset_avg.FonF = mean([offset.FonF]);
offset_avg.FonB = mean([offset.FonB]);
offset_avg.BonF = mean([offset.BonF]);
offset_avg.BonB = mean([offset.BonB]);

% draw sigmoids with the average parameters
xspan = -50:0.1:50;
sigm_avg.FonF = sigmf(xspan,[slope_avg.FonF offset.FonF]);
sigm_avg.FonB = sigmf(xspan,[slope_avg.FonB offset.FonB]);
sigm_avg.BonF = sigmf(xspan,[slope_avg.BonF offset.BonF]);
sigm_avg.BonB = sigmf(xspan,[slope_avg.BonB offset.BonB]);

sigm_avg.congru = sigmf(xspan,[mean([slope_avg.FonF, slope_avg.BonB]), mean([offset.FonF, offset.BonB])]);
sigm_avg.incong = sigmf(xspan,[mean([slope_avg.FonB, slope_avg.BonF]), mean([offset.FonB, offset.BonF])]);

figure('color',[1 1 1])
subplot(4,4,[1 2 5 6])
plot(xspan,sigm_avg.FonF,'linewidth',5,'color',hist_fire_color), hold on % based on value assigned during fire
set(gca,'fontsize',21,'ytick',[],'xtick',[]),ylabel('P(Choose right)'),ylim([-0.025 1.025])

subplot(4,4,[3 4 7 8])
plot(xspan,sigm_avg.BonF,'linewidth',5,'color',hist_boat_color), hold on % based on value assigned during fire
set(gca,'fontsize',21,'ytick',[],'xtick',[]),ylim([-0.025 1.025])

subplot(4,4,[9 10 13 14])
plot(xspan,sigm_avg.FonB,'linewidth',5,'color',hist_fire_color), hold on % based on value assigned during fire
set(gca,'fontsize',21,'ytick',[],'xtick',[]),ylabel('P(Choose right)'),xlabel('v_{right} - v_{left}'),ylim([-0.025 1.025])

subplot(4,4,[11 12 15 16])
plot(xspan,sigm_avg.BonB,'linewidth',5,'color',hist_boat_color), hold on % based on value assigned during fire
set(gca,'fontsize',21,'ytick',[],'xtick',[]),xlabel('v_{right} - v_{left}'),ylim([-0.025 1.025])

% draw sigmoids congruent/incongruent
figure('color',[1 1 1])
subplot(1,2,1)
plot(xspan,sigm_avg.congru,'linewidth',5,'color',hist_fire_color)
set(gca,'fontsize',21,'ytick',[],'xtick',[]),ylim([-0.025 1.025])

subplot(1,2,2)
plot(xspan,sigm_avg.incong,'linewidth',5,'color',hist_boat_color)
set(gca,'fontsize',21,'ytick',[],'xtick',[]),ylim([-0.025 1.025])

clear xspan slope_avg offset
%% separate worst subjects from best subjects
subjectPerformance = mean([slopes.FonF;slopes.BonB],1);
keyboard

%% lme

% create table
tblChoice = table(allSubs,allCond,alldV,alldV_opp,alldF,allCho);
tblChoice.Properties.VariableNames = {'sub','cond','dV','dVo','dF','choice'};

lmeChoice = fitlme(tblChoice,'choice ~ 1 + dV + dVo +  cond + dF + (1 + dV + dVo + cond + dF | sub)');

% plot coefficients
k_cond = lmeChoice.Coefficients(2,2).Estimate;
k_dV   = lmeChoice.Coefficients(3,2).Estimate;
k_dVo  = lmeChoice.Coefficients(4,2).Estimate;
k_dF   = lmeChoice.Coefficients(5,2).Estimate;
ce_cond = lmeChoice.Coefficients(2,3).SE;
ce_dV  = lmeChoice.Coefficients(3,3).SE;
ce_dVo = lmeChoice.Coefficients(4,3).SE;
ce_dF  = lmeChoice.Coefficients(5,3).SE;

figure('color',[1 1 1]),hold on
xspan = -10:0.1:10;
plot(xspan,zeros(length(xspan)),'linestyle',':','color',[0.5 0.5 0.5],'linewidth',2)
set(gca,'fontsize',16,'xtick',[1 2 3 4],'xticklabel',{'Goal','dV same','dV ~same','dF'}),xlim([0.5 4.5]),ylim([-1.5 2])
ylabel('Coefficient')

% errorbar([1 2 3 4],[k_cond k_dV k_dVo k_dF],[ce_cond ce_dV ce_dVo ce_dF],'linestyle','none',...
%     'color','k','linewidth',1.5,'capsize',0),hold on
scatter([1 2 3 4],[k_cond k_dV k_dVo k_dF],180,'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0 0 0])


