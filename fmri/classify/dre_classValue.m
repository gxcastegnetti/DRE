%% dre_classValue
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'class_v0';

% ROIs here are:
%   - ACC
%   - insulae
%   - pHPC
%   - lOFC

%% directories
fs         = filesep;
dir.geoCod = pwd;
idcs       = strfind(dir.geoCod,'/');
dir.dre    = dir.geoCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsaCod = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'class'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];

% paths
addpath([pwd,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

% subs = [4 5 7];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

%% load response patterns and apply mask
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23'};

roiNames = {'sphere_10--20_-54_-8'};
% roiNames = {'box_w-16_16_16-0_36_23'};
% roiNames = {'l_hpc'};
roiNames = {'ofc'};

% apply two masks: one for grey matter, one for ROI
for r = 1:length(roiNames)
    for s = 1:length(subs)
        subjName = ['SF',num2str(subs(s),'%03d')];
        roiMaskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        gmMaskFile =  [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
        roiMask = spm_read_vols(spm_vol(roiMaskFile));
        gmMask = spm_read_vols(spm_vol(gmMaskFile));
        
        % vectorise it
        roiMask = reshape(roiMask, 1, []);
        gmMask = reshape(gmMask, 1, []);
        
        respPatt_foo = responsePatterns.(subjName)(logical(roiMask),:);
        respPatt.(['roi',num2str(r)]).(subjName) = respPatt_foo(~isnan(respPatt_foo(:,1)),:);
    end
end

roiNames = fieldnames(respPatt);
subNames = fieldnames(respPatt.roi1);
figure
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        disp(['sub#',num2str(subs(s))])
        
        % extract presentation indices
        if taskOrd(s) == 1
            objIdx_F = [bData(subs(s)).imagination(1).objIdx'; bData(subs(s)).imagination(3).objIdx'];
            objIdx_B = [bData(subs(s)).imagination(2).objIdx'; bData(subs(s)).imagination(4).objIdx'];
            objVal_F = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(3).val];
            objVal_B = [bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(4).val];
        elseif taskOrd(s) == 2
            objIdx_F = [bData(subs(s)).imagination(2).objIdx'; bData(subs(s)).imagination(4).objIdx'];
            objIdx_B = [bData(subs(s)).imagination(1).objIdx'; bData(subs(s)).imagination(3).objIdx'];
            objVal_F = [bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(4).val];
            objVal_B = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(3).val];
        end
        
        % sort indices
        [~, objIdx_sort_F] = sort(objIdx_F);
        [~, objIdx_sort_B] = sort(objIdx_B);
        
        % sort values accordingly
        objVal_F = objVal_F(objIdx_sort_F);
        objVal_B = objVal_B(objIdx_sort_B);
        
        % take activation patterns
        X_F = respPatt.(roiNames{r}).(subNames{s})(:,1:120)';
        X_B = respPatt.(roiNames{r}).(subNames{s})(:,121:240)';
        
        % fix nans
%         objVal_F(isnan(objVal_F)) = 50*rand;
%         objVal_B(isnan(objVal_B)) = 50*rand;
        
        % add a constant for univoque determination of median
        Y_F = objVal_F + (0.00000001*(1:120))';
        Y_B = objVal_B + (0.00000001*(1:120))';
        
        % find percentiles
        pl_F_33 = prctile(Y_F,100/3);
        pl_F_66 = prctile(Y_F,200/3);
        pl_B_33 = prctile(Y_B,100/3);
        pl_B_66 = prctile(Y_B,200/3);
        
        X_F_red = X_F(Y_F < pl_F_33 | Y_F > pl_F_66,:);
        Y_F_red = Y_F(Y_F < pl_F_33 | Y_F > pl_F_66);
        Y_F_logic = Y_F_red > pl_F_66;
        
        X_B_red = X_B(Y_B < pl_B_33 | Y_B > pl_B_66,:);
        Y_B_red = Y_B(Y_B < pl_B_33 | Y_B > pl_B_66);
        Y_B_logic = Y_B_red > pl_B_66;
        
        % separate into two levels
%         Y_F_logic = Y_F > median(Y_F);
%         Y_B_logic = Y_B > median(Y_B);
        
        nTrials = 80;
        
        %% svm
        
        % FF
        mdl_FF = fitcsvm(X_F_red,Y_F_logic,'BoxConstraint',1000);
        cvmdl_FF = crossval(mdl_FF);
        
        % BB
        mdl_BB = fitcsvm(X_B_red,Y_B_logic,'BoxConstraint',1000);
        cvmdl_BB = crossval(mdl_BB);
        
        % FB
        %         mdl_FB = fitcsvm(X_F_train,Y_F_train);
        %         foo = sum(mdl_FB.Beta.*X_B_test(:)) + mdl_FB.Bias;
        %         prediction_FB(testTrl) = round(1/(1+exp(-foo)));
        %
        %         % BF
        %         mdl_BF = fitcsvm(X_B_train,Y_B_train);
        %         foo = sum(mdl_BF.Beta.*X_F_test(:)) + mdl_BF.Bias;
        %         prediction_BF(testTrl) = round(1/(1+exp(-foo)));
        
        
        % compute performance
        performance_FF(s) = 1 - kfoldLoss(cvmdl_FF);
        performance_BB(s) = 1 - kfoldLoss(cvmdl_BB);
        %         performance_FB(s) = 1 - kfoldLoss(cvmdl_FB);
        %         performance_BF(s) = 1 - kfoldLoss(cvmdl_BF);
        
        subplot(6,6,s),bar([performance_FF(s),performance_BB(s)]),ylim([0.4 0.6])
        
    end
    
    %     performance_within = mean([performance_FF;performance_BB]);
    %     performance_betwee = mean([performance_FB;performance_BF]);
    
    %     figure('color',[1 1 1])
    %     bar([mean(performance_within),mean(performance_betwee)])
    %     ylim([0.49 0.53])
end