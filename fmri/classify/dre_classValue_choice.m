%% dre_classValue_choice
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'classChoice_v0';

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

% roiNames = {'sphere_10--20_-54_-8'};
roiNames = {'sphere_9--28_34_-19'};
roiNamesTrue = roiNames;
roiNames = {'lingual'};
roiNames = {'l_hpc'};

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
        
        clear subjName roiMaskFile gmMaskFile roiMask gmMask respPatt_foo
    end
end, clear r s

roiNames = fieldnames(respPatt);
subNames = fieldnames(respPatt.roi1);

for r = 1:length(roiNames)
    h = figure;
    for s = 1:length(subs)
        
        disp(['sub#',num2str(subs(s))])
        
        % extract presentation indices
        if taskOrd(s) == 1
            objVal_F = [bData(subs(s)).choice(1).chMunc;  bData(subs(s)).choice(3).chMunc];
            objVal_B = [bData(subs(s)).choice(2).chMunc;  bData(subs(s)).choice(4).chMunc];
        elseif taskOrd(s) == 2
            objVal_F = [bData(subs(s)).choice(2).chMunc;  bData(subs(s)).choice(4).chMunc];
            objVal_B = [bData(subs(s)).choice(1).chMunc;  bData(subs(s)).choice(3).chMunc];
        end
        
        % take activation patterns
        X_F = respPatt.(roiNames{r}).(subNames{s})(:,1:48)';
        X_B = respPatt.(roiNames{r}).(subNames{s})(:,49:96)';
        
        % add a constant for univoque determination of median
        Y_F = objVal_F + (0.00000001*(1:48))';
        Y_B = objVal_B + (0.00000001*(1:48))';

        Y_F_logic = Y_F > median(Y_F);
        Y_B_logic = Y_B > median(Y_B);
        
        nTrials = length(Y_F_logic);
        
        clear objVal_F objVal_B
        
        %% optimise models
        Mdl_F = fitcsvm(X_F,Y_F_logic,'OptimizeHyperparameters','auto',...
            'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
            'expected-improvement-plus','ShowPlots',false,'UseParallel',true));
        bc_F = Mdl_F.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
        ks_F = Mdl_F.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
        
        Mdl_B = fitcsvm(X_B,Y_B_logic,'OptimizeHyperparameters','auto',...
            'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
            'expected-improvement-plus','ShowPlots',false,'UseParallel',true));
        bc_B = Mdl_B.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
        ks_B = Mdl_B.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
        
        clear Mdl_F Mdl_B
        
        %% CV
        nSweeps = 100;
        for k = 1:nSweeps
            
            c_F = cvpartition(Y_F_logic,'holdOut',0.1);
            c_B = cvpartition(Y_B_logic,'holdOut',0.1);
            
            idxTrain_F = training(c_F);
            idxTrain_B = training(c_B);
            
            idxTest_F = ~idxTrain_F;
            idxTest_B = ~idxTrain_B;
            
            % train F
            XTrain_F = X_F(idxTrain_F,:);
            yTrain_F = Y_F_logic(idxTrain_F);
            
            % train B
            XTrain_B = X_B(idxTrain_B,:);
            yTrain_B = Y_B_logic(idxTrain_B);
            
            % test F and B
            XTest_same_F = X_F(idxTest_F,:);
            yTest_same_F = Y_F_logic(idxTest_F);
            XTest_same_B = X_B(idxTest_B,:);
            yTest_same_B = Y_B_logic(idxTest_B);
            XTest_diff_F = X_F(idxTest_B,:);
            yTest_diff_F = Y_F_logic(idxTest_B);
            XTest_diff_B = X_B(idxTest_F,:);
            yTest_diff_B = Y_B_logic(idxTest_F);
            
            mdl_F = fitcsvm(XTrain_F,yTrain_F,'BoxConstraint',bc_F,'KernelScale',ks_F);
            mdl_B = fitcsvm(XTrain_B,yTrain_B,'BoxConstraint',bc_B,'KernelScale',ks_B);
            
            % test within and between goals
            label_FF = predict(mdl_F,XTest_same_F);
            label_BB = predict(mdl_B,XTest_same_B);
            label_FB = predict(mdl_F,XTest_diff_B);
            label_BF = predict(mdl_B,XTest_diff_F);
            
            acc_foo_FF(k) = mean(label_FF(:) == yTest_same_F(:));
            acc_foo_BB(k) = mean(label_BB(:) == yTest_same_B(:));
            acc_foo_FB(k) = mean(label_FB(:) == yTest_diff_B(:));
            acc_foo_BF(k) = mean(label_BF(:) == yTest_diff_F(:));
            
        end
        
        % compute performance
        acc_FF(s) = mean(acc_foo_FF);
        acc_BB(s) = mean(acc_foo_BB);
        acc_FB(s) = mean(acc_foo_FB);
        acc_BF(s) = mean(acc_foo_BF);
        
        figure(h)
        subplot(6,6,s)
        bar([1,2],[acc_FF(s),acc_BB(s)],'facecolor',[0.15 0.45 0.75]),hold on
        bar([3.5,4.5],[acc_FB(s),acc_BF(s)],'facecolor',[0.55 0.55 0.55])
        set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',11)
        plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
        ylim([0.4 0.6]),xlim([0 5.5])
        
        clear acc_foo_FF acc_foo_BB acc_foo_FB acc_foo_BF label_FF label_BB label_FB label_BF c_F c_B
        clear XTest_F XTest_B XTrain_F XTrain_B bc_F bc_B ks_F ks_B idxTest_F idxTest_B idxTrain_F idxTrain_B
    end
    
    acc_FF_mean(r) = mean(acc_FF);
    acc_BB_mean(r) = mean(acc_BB);
    acc_FB_mean(r) = mean(acc_FB);
    acc_BF_mean(r) = mean(acc_BF);
    
    figure('color',[1 1 1])
    bar([1,2],[acc_FF_mean(r),acc_BB_mean(r)],'facecolor',[0.15 0.45 0.75]),hold on
    bar([3.5,4.5],[acc_FB_mean(r),acc_BF_mean(r)],'facecolor',[0.55 0.55 0.55])
    set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',14)
    title(roiNamesTrue{r},'fontsize',18)
    plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
    ylim([0.4 0.6]),xlim([0 5.5])
    
end, clear r k s

clear responsePatterns

aaa=mean([acc_FF;acc_BB]);
[h,p,ci,stats] = ttest(aaa-0.5)