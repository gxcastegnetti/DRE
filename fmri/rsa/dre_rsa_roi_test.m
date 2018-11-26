%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_path';

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
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

% roi names
% roiNames = {'sphere_6--4_-52_13','sphere_6--4_-51_20','sphere_6--4_-47_26','sphere_6--4_-41_31','sphere_6--4_-35_34',...
%     'sphere_6--4_-28_34','sphere_6--4_-20_37','sphere_6--4_-13_37','sphere_6--4_-6_39','sphere_6--4_2_38',...
%     'sphere_6--4_9_37','sphere_6--4_13_34','sphere_6--4_20_31','sphere_6--4_27_26','sphere_6--4_34_22',...
%     'sphere_6--4_37_16','sphere_6--4_40_9'};

% roiNames = {'sphere_8--4_-57_15','sphere_8--4_-52_29','sphere_8--4_-39_35','sphere_8--4_-25_39','sphere_8--4_-11_41',...
%     'sphere_8--4_3_40','sphere_8--4_17_34','sphere_8--4_29_26','sphere_8--4_40_15'};

% roiNames = {'sphere_6--4_-40_48','sphere_6--4_-28_46','sphere_6--4_-16_46','sphere_6--4_-4_42','sphere_6--4_8_38',...
%     'sphere_6--4_20_31','sphere_6--4_32_24','sphere_6--4_44_12'};

% roiNames = {'sphere_9-0_-61_25','sphere_9-0_-43_35','sphere_9-0_-25_39','sphere_9-0_-7_40',...
%     'sphere_9-0_11_37','sphere_9-0_28_29','sphere_9-0_41_20','sphere_9-0_47_11'};
% 
% roiNames = {'box_w-16_16_16-0_-61_25','box_w-16_16_16-0_-43_35','box_w-16_16_16-0_-25_39','box_w-16_16_16-0_-7_40',...
%     'box_w-16_16_16-0_11_37','box_w-16_16_16-0_28_29','box_w-16_16_16-0_41_20','box_w-16_16_16-0_47_11'};

roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23'};

roiNames = {'lingual'};

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

%% construct RDMs
RDMs_data = constructRDMs(respPatt, 'SPM', userOptions);
RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');

%% plot RDMs
figureRDMs(RDM_average,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];
for s = 1:length(subs)
    RDMs_model(1,s).name = 'val';
    RDMs_model(1,s).RDM = RDMs_models{s}.val; %#ok<*SAGROW>
    RDMs_model(1,s).color = [0 1 0];
    RDMs_model(2,s).name = 'con';
    RDMs_model(2,s).RDM = RDMs_models{s}.con;
    RDMs_model(2,s).color = [0 1 0];
    RDMs_model(3,s).name = 'fam';
    RDMs_model(3,s).RDM = RDMs_models{s}.fam;
    RDMs_model(3,s).color = [0 1 0];
    RDMs_model(4,s).name = 'oid';
    RDMs_model(4,s).RDM = 1-mat_ID;
    RDMs_model(4,s).color = [0 1 0];
    RDMs_model(5,s).name = 'cxt';
    RDMs_model(5,s).RDM = RDMs_models{s}.cxt;
    RDMs_model(5,s).color = [0 1 0];
end

%% for every region and sub, correlate RDM and model
scoreNames = {'val','con','fam','ID','cxt'};
% h = figure;
for r = 1:length(roiNames)
    for s = 1:size(RDMs_data,2)
        for m = 1:size(RDMs_model,1)
            a = vectorizeRDM(RDMs_data(r,s).RDM);
            b = vectorizeRDM(RDMs_model(m,s).RDM);
            rL2(m,s,r) = corr(b',a','rows','pairwise','type','Spearman');
            rL2(m,s,r) = fisherTransform(rL2(m,s,r));
        end
    end
%     figure(h),subplot(4,5,r)
%     bar(mean(rL2(:,:,r),2)),set(gca,'xticklabel',scoreNames)
%     set(gca,'fontsize',16)
%     ylim([-0.01 0.01]),title(roiNames{r})
end

sphe.val = mean(squeeze(rL2(1,:,:)));
sphe.con = mean(squeeze(rL2(2,:,:)));
sphe.fam = mean(squeeze(rL2(3,:,:)));
sphe.oid = mean(squeeze(rL2(4,:,:)));
sphe.cxt = mean(squeeze(rL2(5,:,:)));

% linear models
mdlVal = fitlm(1:length(sphe.val),sphe.val);
mdlCon = fitlm(1:length(sphe.con),sphe.con);
mdlFam = fitlm(1:length(sphe.fam),sphe.fam);
mdlOid = fitlm(1:length(sphe.oid),sphe.oid);
mdlCxt = fitlm(1:length(sphe.cxt),sphe.cxt);

%% plot linear fits
xspan = 1:length(sphe.val);
figure('color',[1 1 1])

subplot(1,3,1),hold on
plot(xspan,sphe.val,'marker','.','markersize',30,'linestyle','none'),title('Value')
xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
plot(xspan,xspan*mdlVal.Coefficients{2,1}+mdlVal.Coefficients{1,1},'linestyle','--','color','k')

subplot(1,3,2),hold on
plot(xspan,sphe.oid,'marker','.','markersize',30,'linestyle','none'),title('Object ID')
xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
plot(xspan,xspan*mdlOid.Coefficients{2,1}+mdlOid.Coefficients{1,1},'linestyle','--','color','k')

subplot(1,3,3),hold on
plot(xspan,sphe.cxt,'marker','.','markersize',30,'linestyle','none'),title('Context')
xlabel('ROI number'),ylabel('r'),set(gca,'fontsize',16,'xtick',1:7,'xticklabel',1:7)
plot(xspan,xspan*mdlCxt.Coefficients{2,1}+mdlCxt.Coefficients{1,1},'linestyle','--','color','k')

keyboard

%% ttest
for r = 1:length(roiNames)
    for m = 1:size(RDMs_model,1)
        scores = rL2(m,:,r);
        [h,p(m,r),~,~] = ttest(scores);
    end
end

%%
keyboard
subsBest = [23 18 5 3 21 11 10 20 17 28 24  1 15 22];
subsWors = [2  14 6 4  7  9 27 26 12 16 19 13  8 25];

bestVsWorsePlot = figure('color',[1 1 1]);
x_err = [0.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15];
for r = 1:size(RDMs_data,1)
    
    % compute mean and sem for best and worse perf subs
    meanBest = mean(rL2(:,subsBest,r),2);
    meanWors = mean(rL2(:,subsWors,r),2);
    semBest = std(rL2(:,subsBest,r)')/sqrt(length(subs));
    semWors = std(rL2(:,subsWors,r)')/sqrt(length(subs));
    
    % plot bars
    figure(bestVsWorsePlot)
    subplot(4,3,r)
    bar([meanBest,meanWors])
    if r == 1,legend('best','worse','location','northwest'),end
    ylim([-0.015 0.009])
    set(gca,'xticklabel',scoreNames,'fontsize',14),hold on
    
    % plot error bars
%     foo_mean = [meanBest';meanWors'];
%     foo_mean = foo_mean(:);
%     foo_err = [semBest',semWors'];
%     foo_err = foo_err(:);
%     errorbar(x_err,foo_mean,foo_err,'linestyle','none','color','k')
end



