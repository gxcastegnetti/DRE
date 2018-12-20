%% dre_rsa_roi
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'test';

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% load response patterns and apply mask
% filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ima/rsaPatterns_sl.mat';
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_choice/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

roiNames = {'l_hpc','r_hpc','l_ofc'};

% apply two masks: one for grey matter, one for ROI
for r = 1:length(roiNames)
    for s = 1:length(subs)
        subjName = ['SF',num2str(subs(s),'%03d')];
        roiMaskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        gmMaskFile = [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
        roiMask = spm_read_vols(spm_vol(roiMaskFile));
        gmMask = spm_read_vols(spm_vol(gmMaskFile));
        
        % vectorise it
        roiMask = reshape(roiMask, 1, []);
        gmMask = reshape(gmMask, 1, []);
        %         totMask = logical(gmMask.*roiMask);
        
        respPatt_foo = responsePatterns.(subjName)(logical(roiMask) & logical(gmMask),:);
        respPatt.(['roi',num2str(r)]).(subjName) = respPatt_foo(~isnan(respPatt_foo(:,1)),:);
    end
end

%% create structures of rearranged response patterns
if size(respPatt_foo,2) == 240
    respPatt_acc2ses = respPatt;
    respPatt_acc2val = respPatt;
    respPatt_acc2fam = respPatt;
    
    for r = 1:length(roiNames)
        for s = 1:length(subs)
            subjName = ['SF',num2str(subs(s),'%03d')];
            respPatt_acc2ses.(['roi',num2str(r)]).(subjName) = respPatt_acc2ses.(['roi',num2str(r)]).(subjName)(:,ordData(subs(s)).norm2sessions);
            respPatt_acc2val.(['roi',num2str(r)]).(subjName) = respPatt_acc2val.(['roi',num2str(r)]).(subjName)(:,ordData(subs(s)).norm2val_disc);
            respPatt_acc2fam.(['roi',num2str(r)]).(subjName) = respPatt_acc2fam.(['roi',num2str(r)]).(subjName)(:,ordData(subs(s)).norm2fam_disc);
        end
    end
end

%% construct RDMs
RDMs_data = constructRDMs(respPatt, 'SPM', userOptions);
RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');
% RDM_average.name = 'Pearson';

%% plot RDMs
figureRDMs(RDM_average,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];

% if imagination
if size(respPatt_foo,2) == 240
    for s = 1:length(subs)
        RDM_model(1,s).name = 'val';
        RDM_model(1,s).RDM = RDMs_models{s}.val; %#ok<*SAGROW>
        RDM_model(1,s).color = [0 1 0];
        RDM_model(2,s).name = 'con';
        RDM_model(2,s).RDM = RDMs_models{s}.con;
        RDM_model(2,s).color = [0 1 0];
        RDM_model(3,s).name = 'fam';
        RDM_model(3,s).RDM = RDMs_models{s}.fam;
        RDM_model(3,s).color = [0 1 0];
        RDM_model(4,s).name = 'oid';
        RDM_model(4,s).RDM = 1-mat_ID;
        RDM_model(4,s).color = [0 1 0];
        RDM_model(5,s).name = 'cxt';
        RDM_model(5,s).RDM = RDMs_models{s}.cxt;
        RDM_model(5,s).color = [0 1 0];
    end
    
    % for every region and sub, correlate RDM and model
    scoreNames = {'Value','Confidence','Familiarity','Object ID','Context'};
    for r = 1:length(roiNames)
        for s = 1:size(RDMs_data,2)
            for m = 1:size(RDMs_data,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
                % rL2(m,s) = fisherTransform(rL2(m,s));
            end
        end
    end
    
    % else, if choice
elseif size(respPatt_foo,2) == 96
    for s = 1:length(subs)
        RDM_model(1,s).name = 'dVal';
        RDM_model(1,s).RDM = RDMs_models{s}.choice.dVal; %#ok<*SAGROW>
        RDM_model(1,s).color = [0 1 0];
        RDM_model(2,s).name = 'cMun';
        RDM_model(2,s).RDM = RDMs_models{s}.choice.cMun;
        RDM_model(2,s).color = [0 1 0];
        RDM_model(3,s).name = 'chos';
        RDM_model(3,s).RDM = RDMs_models{s}.choice.Chos;
        RDM_model(3,s).color = [0 1 0];
        RDM_model(4,s).name = 'unch';
        RDM_model(4,s).RDM = RDMs_models{s}.choice.Unch;
        RDM_model(4,s).color = [0 1 0];
        RDM_model(5,s).name = 'ccxt';
        RDM_model(5,s).RDM = RDMs_models{s}.choice.ccxt;
        RDM_model(5,s).color = [0 1 0];
    end
    
    % for every region and sub, correlate RDM and model
    scoreNames = {'Value diff.','Chosen - unch.','Chosen value','Unch. value','Context'};
    for r = 1:length(roiNames)
        for s = 1:size(RDMs_data,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDMs_data(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
                % rL2(m,s) = fisherTransform(rL2(m,s));
            end
        end
    end
end

%% interpolate

% sphe.val = mean(squeeze(rL2(1,:,:)));
% sphe.con = mean(squeeze(rL2(2,:,:)));
% sphe.fam = mean(squeeze(rL2(3,:,:)));
% sphe.oid = mean(squeeze(rL2(4,:,:)));
% sphe.cxt = mean(squeeze(rL2(5,:,:)));
%
% % linear models
% mdlVal = fitlm(1:length(sphe.val),sphe.val);
% mdlCon = fitlm(1:length(sphe.con),sphe.con);
% mdlFam = fitlm(1:length(sphe.fam),sphe.fam);
% mdlOid = fitlm(1:length(sphe.oid),sphe.oid);
% mdlCxt = fitlm(1:length(sphe.cxt),sphe.cxt);
%
% %% plot linear fits
% xspan = 1:length(sphe.val);
% figure('color',[1 1 1])
%
% subplot(1,3,1),hold on
% plot(xspan,sphe.val,'marker','.','markersize',30,'linestyle','none'),title('Value')
% xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
% plot(xspan,xspan*mdlVal.Coefficients{2,1}+mdlVal.Coefficients{1,1},'linestyle','--','color','k')
%
% subplot(1,3,2),hold on
% plot(xspan,sphe.oid,'marker','.','markersize',30,'linestyle','none'),title('Object ID')
% xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
% plot(xspan,xspan*mdlOid.Coefficients{2,1}+mdlOid.Coefficients{1,1},'linestyle','--','color','k')
%
% subplot(1,3,3),hold on
% plot(xspan,sphe.cxt,'marker','.','markersize',30,'linestyle','none'),title('Context')
% xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
% plot(xspan,xspan*mdlCxt.Coefficients{2,1}+mdlCxt.Coefficients{1,1},'linestyle','--','color','k')
%
% keyboard

%% ttest
for r = 1:length(roiNames)
    for m = 1:size(RDM_model,1)
        scores = corrRoiModel(r,:,m);
        [h,p(m,r),~,~] = ttest(scores);
    end
end

%% mean and sem
means = squeeze(mean(corrRoiModel,2));
sems  = squeeze(std(corrRoiModel,0,2)/sqrt(numel(subs)));

%% plot
figure('color',[1 1 1])
hb = bar(means); hold on
% myColors = [0,0,1;1,0,0];
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData',means(:,ib),sems(:,ib),'k.')
    %     hb(ib).FaceColor = myColors(ib,:);
end

legend(scoreNames,'location','northwest'),set(gca,'fontsize',18,'xtick',1:4,...
    'xticklabels',{'HPC','Cingul.','vmPFC','OFC'})
ylabel('Correlation(ROI, model)')

%%
% keyboard
% subsBest = [23 18 5 3 21 11 10 20 17 28 24  1 15 22];
% subsWors = [2  14 6 4  7  9 27 26 12 16 19 13  8 25];
%
% bestVsWorsePlot = figure('color',[1 1 1]);
% x_err = [0.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15];
% for r = 1:size(RDMs_data,1)
%
%     % compute mean and sem for best and worse perf subs
%     meanBest = mean(rL2(r,subsBest,r),2);
%     meanWors = mean(rL2(:,subsWors,r),2);
%     semBest = std(rL2(:,subsBest,r)')/sqrt(length(subs));
%     semWors = std(rL2(:,subsWors,r)')/sqrt(length(subs));
%
%     % plot bars
%     figure(bestVsWorsePlot)
%     subplot(4,3,r)
%     bar([meanBest,meanWors])
%     if r == 1,legend('best','worse','location','northwest'),end
%     ylim([-0.015 0.009])
%     set(gca,'xticklabel',scoreNames,'fontsize',14),hold on
%
%     % plot error bars
% %     foo_mean = [meanBest';meanWors'];
% %     foo_mean = foo_mean(:);
% %     foo_err = [semBest',semWors'];
% %     foo_err = foo_err(:);
% %     errorbar(x_err,foo_mean,foo_err,'linestyle','none','color','k')
% end



