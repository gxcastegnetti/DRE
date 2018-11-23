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
% subs = subs(6);
% taskOrd = taskOrd(6);
% subs = [4 5 7 8 9 13 14];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

%% load response patterns and apply mask
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23'};

% roiNames = {'sphere_10--20_-54_-8'};
% roiNames = {'box_w-16_16_16-0_20_36'};
% roiNames = {'l_hpc'};
roiNames = {'sphere_9--28_34_-19'};

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
        
        respPatt_foo = responsePatterns.(subjName)(logical(roiMask) & logical(gmMask),:);
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
        
        for i_F = 1:size(X_F)
            X_F(i_F,:) = zscore(X_F(i_F,:));
            X_B(i_F,:) = zscore(X_B(i_F,:));
        end
        
        % fix nans
        objVal_F(isnan(objVal_F)) = 50*rand;
        objVal_B(isnan(objVal_B)) = 50*rand;
        
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
%         [coef, info] = lassoglm(X_F_red,Y_F_logic,'binomial','CV',10,'LambdaRatio',0.2,'NumLambda',50);
%         lassoPlot(coef,info,'plottype','CV');
%         legend('show') % Show legend
%         
        nSweeps = 25;
        acc_FF_foo = nan(nSweeps,1);
        holdOutFraction = 0.1;
        for k = 1:nSweeps
            
            c_F = cvpartition(Y_F_logic,'HoldOut',holdOutFraction);
            c_B = cvpartition(Y_B_logic,'HoldOut',holdOutFraction);
            
            idxTrain_F = training(c_F);
            idxTrain_B = training(c_B);
            
            idxTest_F = ~idxTrain_F;
            idxTest_B = ~idxTrain_B;
            
            XTrain_F = X_F_red(idxTrain_F,:);
            yTrain_F = Y_F_logic(idxTrain_F);
            XTrain_B = X_B_red(idxTrain_B,:);
            yTrain_B = Y_B_logic(idxTrain_B);
            
            XTest_F = X_F_red(idxTest_F,:);
            yTest_F = Y_F_logic(idxTest_F);
            XTest_B = X_B_red(idxTest_B,:);
            yTest_B = Y_B_logic(idxTest_B);
            
            [b_F,i_F] = lassoglm(XTrain_F,yTrain_F,'binomial','lambda',0.01);
            [b_B,i_B] = lassoglm(XTrain_F,yTrain_F,'binomial','lambda',0.01);
            
            for j = 1:round(holdOutFraction*nTrials)
                
                % take test samples
                XTest_F_foo = XTest_F(j,:);
                XTest_B_foo = XTest_B(j,:);
                
                % FF
                foo = sum(b_F.*XTest_F_foo(:)) + i_F.Intercept;
                label_FF(j) = round(1./(1+exp(-foo)));
                
                % BB
                foo = sum(b_B.*XTest_B_foo(:)) + i_B.Intercept;
                label_BB(j) = round(1./(1+exp(-foo)));
                
                % FB
                foo = sum(b_F.*XTest_B_foo(:)) + i_F.Intercept;
                label_FB(j) = round(1./(1+exp(-foo)));
                
                % BF
                foo = sum(b_B.*XTest_F_foo(:)) + i_B.Intercept;
                label_BF(j) = round(1./(1+exp(-foo)));
            end
            
            acc_FF_foo(k) = mean(label_FF(:) == yTest_F);
            acc_BB_foo(k) = mean(label_BB(:) == yTest_B);
            acc_FB_foo(k) = mean(label_FB(:) == yTest_F);
            acc_BF_foo(k) = mean(label_BF(:) == yTest_B);
        end
        
        acc_FF_mean(s) = mean(acc_FF_foo);
        acc_BB_mean(s) = mean(acc_BB_foo);
        acc_FB_mean(s) = mean(acc_FB_foo);
        acc_BF_mean(s) = mean(acc_BF_foo);
        
        %         clear prediction_FF prediction_BB objIdx_F objIdx_B
        
    end
    figure('color',[1 1 1]),bar([mean(acc_FF_mean),mean(acc_BB_mean);mean(acc_FB_mean),mean(acc_BF_mean)])
    ylim([0.45 0.55])
    aaa = acc_FF_mean + acc_BB_mean - acc_FB_mean - acc_BF_mean;
end