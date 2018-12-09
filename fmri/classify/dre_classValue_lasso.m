%% dre_classValue
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'class_v0';
betaid       = 'rsa_pulse_ima';

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
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ima/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23'};

% roiNames = {'sphere_10--20_-54_-8'};
% roiNames = {'box_w-16_16_16-0_20_36'};
% roiNames = {'l_hpc'};
roiNames = {'midOcc','lingual','imaginationValue','lp_hpc','rp_hpc','la_hpc','ra_hpc','pcc','mcc','acc','ofc'};
roiNames = {'ofc'};
% roiNames = {'lingual'};
roiNamesTrue = roiNames;

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
        
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        % take only those corresponding to conditions
        B = real(B([1:60,67:126,133:192,199:258],:));
        B = B(toNormalOrder,:)';
        
        respPatt.(['roi',num2str(r)]).(subjName) = B;

    end
end

roiNames = fieldnames(respPatt);
subNames = fieldnames(respPatt.roi1);
hMean = figure('color',[1 1 1]);
hSub = figure('color',[1 1 1]);
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
        Y_F = objVal_F + (0.00001*(1:120))';
        Y_B = objVal_B + (0.00001*(1:120))';
        
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
        
        %% logistic regression
        
        % FF
        %         [~, info_F] = lassoglm(X_F_red,Y_F_logic,'binomial','CV',10,'LambdaRatio',0.2,'NumLambda',50);
        %         lassoPlot(coef,info,'plottype','CV');
        %         legend('show') % Show legend
        
        % FF
        %         [~, info_B] = lassoglm(X_B_red,Y_B_logic,'binomial','CV',10,'LambdaRatio',0.2,'NumLambda',50);
        %         lassoPlot(coef,info,'plottype','CV');
        %         legend('show') % Show legend
        
        numFolds = 10;
        for k = 1:numFolds
            
            c_F = cvpartition(Y_F_logic,'kFold',numFolds);
            c_B = cvpartition(Y_B_logic,'kFold',numFolds);
            
            idxTrain_F = training(c_F,k);
            idxTrain_B = training(c_B,k);
            
            idxTest_F = ~idxTrain_F;
            idxTest_B = ~idxTrain_B;
            
            XTrain_F = X_F_red(idxTrain_F,:);
            yTrain_F = Y_F_logic(idxTrain_F);
            XTrain_B = X_B_red(idxTrain_B,:);
            yTrain_B = Y_B_logic(idxTrain_B);
            
            XTest_same_F = X_F_red(idxTest_F,:);
            yTest_same_F = Y_F_logic(idxTest_F);
            XTest_same_B = X_B_red(idxTest_B,:);
            yTest_same_B = Y_B_logic(idxTest_B);
            XTest_diff_F = X_F_red(idxTest_B,:);
            yTest_diff_F = Y_F_logic(idxTest_B);
            XTest_diff_B = X_B_red(idxTest_F,:);
            yTest_diff_B = Y_B_logic(idxTest_F);
            
            [b_F,i_F] = lassoglm(XTrain_F,yTrain_F,'binomial','lambda',0.001);
            [b_B,i_B] = lassoglm(XTrain_B,yTrain_B,'binomial','lambda',0.001);
            %             [b_F,i_F] = lassoglm(XTrain_F,yTrain_F,'binomial','lambda',info_F.LambdaMinDeviance);
            %             [b_B,i_B] = lassoglm(XTrain_B,yTrain_B,'binomial','lambda',info_B.LambdaMinDeviance);
            
            for j = 1:round(nTrials/numFolds)
                
                % take test samples
                XTest_same_F_foo = XTest_same_F(j,:);
                XTest_same_B_foo = XTest_same_B(j,:);
                XTest_diff_B_foo = XTest_diff_B(j,:);
                XTest_diff_F_foo = XTest_diff_F(j,:);
                
                % FF
                foo = sum(b_F.*XTest_same_F_foo(:)) + i_F.Intercept;
                label_FF(j) = round(1./(1+exp(-foo)));
                
                % BB
                foo = sum(b_B.*XTest_same_B_foo(:)) + i_B.Intercept;
                label_BB(j) = round(1./(1+exp(-foo)));
                
                % FB
                foo = sum(b_F.*XTest_diff_B_foo(:)) + i_F.Intercept;
                label_FB(j) = round(1./(1+exp(-foo)));
                
                % BF
                foo = sum(b_B.*XTest_diff_F_foo(:)) + i_B.Intercept;
                label_BF(j) = round(1./(1+exp(-foo)));
            end
            
            acc_foo_FF(k) = mean(label_FF(:) == yTest_same_F(:));
            acc_foo_BB(k) = mean(label_BB(:) == yTest_same_B(:));
            acc_foo_FB(k) = mean(label_FB(:) == yTest_diff_B(:));
            acc_foo_BF(k) = mean(label_BF(:) == yTest_diff_F(:));
        end
        
        acc_FF_meanSub(s) = mean(acc_foo_FF);
        acc_BB_meanSub(s) = mean(acc_foo_BB);
        acc_FB_meanSub(s) = mean(acc_foo_FB);
        acc_BF_meanSub(s) = mean(acc_foo_BF);
        
        %         clear prediction_FF prediction_BB objIdx_F objIdx_B
        
        figure(hSub)
        subplot(6,6,s)
        bar([1,2],[acc_FF_meanSub(s),acc_BB_meanSub(s)],'facecolor',[0.15 0.45 0.75]),hold on
        bar([3.5,4.5],[acc_FB_meanSub(s),acc_BF_meanSub(s)],'facecolor',[0.55 0.55 0.55])
        set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',11)
        plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
        ylim([0.4 0.6]),xlim([0 5.5])
        
        
    end
    acc_FF_meanRoi(r) = mean(acc_FF_meanSub);
    acc_BB_meanRoi(r) = mean(acc_BB_meanSub);
    acc_FB_meanRoi(r) = mean(acc_FB_meanSub);
    acc_BF_meanRoi(r) = mean(acc_BF_meanSub);
    
    figure(hMean),subplot(2,3,r)
    bar([1,2],[acc_FF_meanRoi(r),acc_BB_meanRoi(r)],'facecolor',[0.15 0.45 0.75]),hold on
    bar([3.5,4.5],[acc_FB_meanRoi(r),acc_BF_meanRoi(r)],'facecolor',[0.55 0.55 0.55])
    set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',14)
    title(roiNamesTrue{r},'fontsize',18)
    plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
    ylim([0.4 0.6]),xlim([0 5.5])
    
    [~,pSame(r),~,~] = ttest(mean([acc_FF_meanSub(:);acc_BB_meanSub(:)],2)-0.5)
    [~,pOppo(r),~,~] = ttest(mean([acc_FB_meanSub(:);acc_BF_meanSub(:)],2)-0.5)
    
end