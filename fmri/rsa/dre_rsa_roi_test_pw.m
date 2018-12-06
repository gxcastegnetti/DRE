%% dre_rsa_roi_test_pw
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_pw';
betaid       = 'rsa_pulse_ima';

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
subs = [4 5 7:9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% define betas and apply mask
roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23'};
% roiNames = {'mask_sl_val_l_lg','mask_sl_val_r_lg','mask_sl_val_rm_lg','mask_sl_val_rma_lg',...
%     'l_hpc','r_hpc','mask_sl_val_lp_ins','mask_sl_val_rp_ins','mask_sl_val_ra_ins','mask_sl_val_l_mcc','mask_sl_val_l_acc','mask_sl_val_l_ofc'};

% roiNames = {'lingual','itc','phpc','hpc','angular','ins_la','pcc','mcc','caudate','putamen','pfc_vm','subgenual','ofc'};

roiNames = {'lingual','lp_hpc','rp_hpc','la_hpc','ra_hpc','pcc','mcc','pfc_vm','ofc'};
roiNames = {'midOcc','lingual','phpc','itc','lp_hpc','rp_hpc','la_hpc','ra_hpc','pcc','mcc','acc','ins_la','ofc'};
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

% apply two masks: one for grey matter, one for ROI
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        if taskOrd(s) == 1
            for obj = 1:120
                objIdx_sess_F = [bData(subs(s)).imagination(1).objIdx'; nan(60,1); bData(subs(s)).imagination(3).objIdx'; nan(60,1)];
                objIdx_sess_B = [nan(60,1); bData(subs(s)).imagination(2).objIdx'; nan(60,1); bData(subs(s)).imagination(4).objIdx'];
                objIdx_F(obj) = find(objIdx_sess_F == obj);
                objIdx_B(obj) = find(objIdx_sess_B == obj);
            end
        elseif taskOrd(s) == 2
            for obj = 1:120
                objIdx_sess_B = [bData(subs(s)).imagination(1).objIdx'; nan(60,1); bData(subs(s)).imagination(3).objIdx'; nan(60,1)];
                objIdx_sess_F = [nan(60,1); bData(subs(s)).imagination(2).objIdx'; nan(60,1); bData(subs(s)).imagination(4).objIdx'];
                objIdx_F(obj) = find(objIdx_sess_F == obj);
                objIdx_B(obj) = find(objIdx_sess_B == obj);
            end
        end
        toNormalOrder = [objIdx_F,objIdx_B];
        
        % SPM file from 1st level analysis
        subjSPMFile = [dir.beta,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        % whiten betas
        B = noiseNormaliseBeta_roi(SPM,subjMaskFile);
        
        % take only those corresponding to conditions
        if size(B,1) == 124
            B = real(B([1:24,31:54,61:84,91:114],:));
        elseif size(B,1) == 268
            B = real(B([1:60,67:126,133:192,199:258],:));
            B = B(toNormalOrder,:);
        end
        
        fooDir = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
        save([fooDir,fs,'toNormalOrder',fs,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
        %         B = B(ordData(subs(s)).norm2val_cont,:);
        
        % construct RDM
        rdm = squareform(pdist(B,'correlation'));
        %     RDM_struct(s).RDM = rdm;
        %     RDM_struct(s).name = ['sub#',num2str(subs(s),'%03d')];
        %     RDM_struct(s).color = [0 0 1];
        RDM_brain(r,s).RDM = rdm;
    end
end

% figureRDMs(RDM_struct,userOptions)

% RDM_mean.RDM = mean(RDM_brain,3);
% RDM_mean.name = 'Mahalanobis';
% RDM_mean.color = [0 0 1];
% figureRDMs(RDM_mean,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];

% if imagination
if size(rdm,1) == 240
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
        for s = 1:size(RDM_brain,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
                % rL2(m,s) = fisherTransform(rL2(m,s));
            end
        end
    end
    
    % else, if choice
elseif size(rdm,1) == 96
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
        for s = 1:size(RDM_brain,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
                % rL2(m,s) = fisherTransform(rL2(m,s));
            end
        end
    end
end

%% ttest
for r = 1:length(roiNames)
    for m = 1:size(RDM_model,1)
        scores = corrRoiModel(r,:,m);
        [h,p(r,m),~,~] = ttest(scores,0,'Tail','right');
    end
end

%% mean and sem
means = squeeze(mean(corrRoiModel,2));
sems  = squeeze(std(corrRoiModel,0,2)/sqrt(numel(subs)));

%% plot
roiNamesTrue = {'midOcc','Lingual','phpc','itc','lpHPC','rpHPC','laHPC','raHPC','PCC','MCC','ACC','aIns','OFC'};
% roiNamesTrue = {'1','2','3','4','5','6','7'};
figure('color',[1 1 1])
hb = bar(means); hold on
% myColors = [0,0,1;1,0,0];
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData',means(:,ib),sems(:,ib),'k.')
    %     hb(ib).FaceColor = myColors(ib,:);
end

% legend(scoreNames,'location','northwest'),set(gca,'fontsize',18,'xtick',1:4,...
%     'xticklabels',{'HPC','Cingul.','vmPFC','OFC'})
legend(scoreNames,'location','northwest'),set(gca,'fontsize',18,'xtick',1:numel(roiNamesTrue),...
    'xticklabels',roiNamesTrue)
ylabel('Correlation(ROI, model)')

%% correlation between representations
for r1 = 1:length(roiNames)
    for r2 = 1:length(roiNames)
        for s = 1:size(RDM_brain,2)
            roi1 = vectorizeRDM(RDM_brain(r1,s).RDM);
            roi2 = vectorizeRDM(RDM_brain(r2,s).RDM);
            corrRoiRoi(r1,r2,s) = corr(roi1(:),roi2(:),'rows','pairwise','type','Spearman');
        end
    end
end
corrRoiRoi_mean = mean(corrRoiRoi,3);
figure,imagesc(corrRoiRoi_mean,[0 0.18])
set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45)

%% Compare good and bad subjects
% subsBest = [25 20 6 4 23 13 11 22 19 30 26  1 17 24];
% subsWors = [2  16 7 5  8 10 29 28 14 18 21 15 10 27];
% 
% corrRoiRoi_mean_best = mean(corrRoiRoi(:,:,subsBest),3);
% figure,imagesc(corrRoiRoi_mean_best,[0 0.18])
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
% 
% corrRoiRoi_mean_wors = mean(corrRoiRoi(:,:,subsWors),3);
% figure,imagesc(corrRoiRoi_mean_wors,[0 0.18])
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
% 
% corrRoiRoi_mean_diff = corrRoiRoi_mean_best - corrRoiRoi_mean_wors;
% figure,imagesc(corrRoiRoi_mean_diff)
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
% 
% for i = 1:numel(roiNames)
%     for j = 1:numel(roiNames)
%         bbb = corrRoiRoi(i,j,subsBest);
%         www = corrRoiRoi(i,j,subsWors);
%         [~,p(i,j),ci,stats] = ttest(bbb-www);
%         ttt = stats.tstat;
%         aaa(i,j) = ttt;
%     end
% end
% figure,imagesc(aaa)
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
