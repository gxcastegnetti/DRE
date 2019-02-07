%% dre_classValue
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'class_v0';
betaid       = 'rsa_pulse_ima';

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
subs = [4 5 7 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40:43 45 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,4)];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

%% load response patterns and apply mask
roiNames = {'rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm'};
roiNamesTrue = roiNames;

%% apply two masks: one for grey matter, one for ROI
% for r = 1:length(roiNames)
%     for s = 1:length(subs)
%         subjName = ['SF',num2str(subs(s),'%03d')];
%         roiMaskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
%         gmMaskFile =  [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
%         roiMask = spm_read_vols(spm_vol(roiMaskFile));
%         gmMask = spm_read_vols(spm_vol(gmMaskFile));
%
%         % vectorise it
%         roiMask = reshape(roiMask, 1, []);
%         gmMask = reshape(gmMask, 1, []);
%
%         respPatt_foo = responsePatterns.(subjName)(logical(roiMask) & logical(gmMask),:);
%         respPatt.(['roi',num2str(r)]).(subjName) = respPatt_foo(~isnan(respPatt_foo(:,1)),:);
%
%         clear subjName roiMaskFile gmMaskFile roiMask gmMask respPatt_foo
%     end
% end, clear r s

%%
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

%%
roiNames = fieldnames(respPatt);
subNames = fieldnames(respPatt.roi1);

figAllRois = figure('color',[1 1 1]);
for r = 1:length(roiNames)
    %     hfig{r} = figure;
    for s = 1:length(subs)
        
        disp(['sub#',num2str(subs(s))])
        
        % extract presentation indices
        if taskOrd(s) == 1
            objIdx_F = [bData(subs(s)).imagination(1).objIdx'; bData(subs(s)).imagination(3).objIdx'];
            objIdx_B = [bData(subs(s)).imagination(2).objIdx'; bData(subs(s)).imagination(4).objIdx'];
            objVal_F = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(3).val];
            objVal_B = [bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(4).val];
            objCon_F = [bData(subs(s)).imagination(1).con; bData(subs(s)).imagination(3).con];
            objCon_B = [bData(subs(s)).imagination(2).con; bData(subs(s)).imagination(4).con];
        elseif taskOrd(s) == 2
            objIdx_F = [bData(subs(s)).imagination(2).objIdx'; bData(subs(s)).imagination(4).objIdx'];
            objIdx_B = [bData(subs(s)).imagination(1).objIdx'; bData(subs(s)).imagination(3).objIdx'];
            objVal_F = [bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(4).val];
            objVal_B = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(3).val];
            objCon_F = [bData(subs(s)).imagination(2).con; bData(subs(s)).imagination(4).con];
            objCon_B = [bData(subs(s)).imagination(1).con; bData(subs(s)).imagination(3).con];
        end
        
        % sort indices
        [~, objIdx_sort_F] = sort(objIdx_F);
        [~, objIdx_sort_B] = sort(objIdx_B);
        
        % sort values and confidences accordingly
        objVal_F = objVal_F(objIdx_sort_F);
        objVal_B = objVal_B(objIdx_sort_B);
        objCon_F = objCon_F(objIdx_sort_F);
        objCon_B = objCon_B(objIdx_sort_B);
        
        % take activation patterns
        X_F = respPatt.(roiNames{r}).(subNames{s})(:,1:120)';
        X_B = respPatt.(roiNames{r}).(subNames{s})(:,121:240)';
        
        % fix nans
        objVal_F(isnan(objVal_F)) = ceil(50*rand);
        objVal_B(isnan(objVal_B)) = ceil(50*rand);
        
        % add a constant for univoque determination of median
        Y_F = objVal_F + (0.00001*(1:120))';
        Y_B = objVal_B + (0.00001*(1:120))';
        %         Y_F = objCon_F + (0.00001*(1:120))';
        %         Y_B = objCon_B + (0.00001*(1:120))';
        
        % find percentiles
%         pl_F_low = prctile(Y_F,100/3);
%         pl_F_hig = prctile(Y_F,200/3);
%         pl_B_low = prctile(Y_B,100/3);
%         pl_B_hig = prctile(Y_B,200/3);
        
        pl_F_low = prctile(Y_F,50);
        pl_F_hig = prctile(Y_F,50);
        pl_B_low = prctile(Y_B,50);
        pl_B_hig = prctile(Y_B,50);
        
        X_F_red = X_F(Y_F < pl_F_low | Y_F > pl_F_hig,:);
        Y_F_red = Y_F(Y_F < pl_F_low | Y_F > pl_F_hig);
        %         X_F_red = X_F;
        %         Y_F_red = Y_F;
        Y_F_logic = Y_F_red > pl_F_hig;
        
        X_B_red = X_B(Y_B < pl_B_low | Y_B > pl_B_hig,:);
        Y_B_red = Y_B(Y_B < pl_B_low | Y_B > pl_B_hig);
        %         X_B_red = X_B;
        %         Y_B_red = Y_B;
        Y_B_logic = Y_B_red > pl_B_hig;
        
        nTrials = 120;
        
        clear pl_F_low pl_B_low pl_F_hig pl_B_hig
        clear objVal_F objVal_B objIdx_sort_F objIdx_sort_B objIdx_F objIdx_B
        
        %% optimise models
        Mdl_F = fitcsvm(X_F_red,Y_F_logic,'OptimizeHyperparameters','auto',...
            'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
            'expected-improvement-plus','ShowPlots',false,'UseParallel',true));
        bc_F = Mdl_F.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
        ks_F = Mdl_F.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
        
        Mdl_B = fitcsvm(X_B_red,Y_B_logic,'OptimizeHyperparameters','auto',...
            'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
            'expected-improvement-plus','ShowPlots',false,'UseParallel',true));
        bc_B = Mdl_B.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;
        ks_B = Mdl_B.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
        
        clear Mdl_F Mdl_B
        
        %% CV
        nSweeps = 1000;
        for k = 1:nSweeps
            
            c_F = cvpartition(Y_F_logic,'holdOut',0.1);
            c_B = cvpartition(Y_B_logic,'holdOut',0.1);
            
            idxTrain_F = training(c_F);
            idxTrain_B = training(c_B);
            
            idxTest_F = ~idxTrain_F;
            idxTest_B = ~idxTrain_B;
            
            % train F
            XTrain_F = X_F_red(idxTrain_F,:);
            yTrain_F = Y_F_logic(idxTrain_F);
            
            % train B
            XTrain_B = X_B_red(idxTrain_B,:);
            yTrain_B = Y_B_logic(idxTrain_B);
            
            % test F and B
            XTest_same_F = X_F_red(idxTest_F,:);
            yTest_same_F = Y_F_logic(idxTest_F);
            XTest_same_B = X_B_red(idxTest_B,:);
            yTest_same_B = Y_B_logic(idxTest_B);
            XTest_diff_F = X_F_red(idxTest_B,:);
            yTest_diff_F = Y_F_logic(idxTest_B);
            XTest_diff_B = X_B_red(idxTest_F,:);
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
        acc_FF(r,s) = mean(acc_foo_FF);
        acc_BB(r,s) = mean(acc_foo_BB);
        acc_FB(r,s) = mean(acc_foo_FB);
        acc_BF(r,s) = mean(acc_foo_BF);
        
        %         figure(hfig{r})
        %         subplot(6,6,s)
        %         bar([1,2],[acc_FF(s),acc_BB(s)],'facecolor',[0.15 0.45 0.75]),hold on
        %         bar([3.5,4.5],[acc_FB(s),acc_BF(s)],'facecolor',[0.55 0.55 0.55])
        %         set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',11)
        %         plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
        %         ylim([0.4 0.6]),xlim([0 5.5])
        
        clear acc_foo_FF acc_foo_BB acc_foo_FB acc_foo_BF label_FF label_BB label_FB label_BF c_F c_B
        clear XTest_F XTest_B XTrain_F XTrain_B bc_F bc_B ks_F ks_B idxTest_F idxTest_B idxTrain_F idxTrain_B
    end
    
    acc_FF_mean(r) = mean(acc_FF(r,:));
    acc_BB_mean(r) = mean(acc_BB(r,:));
    acc_FB_mean(r) = mean(acc_FB(r,:));
    acc_BF_mean(r) = mean(acc_BF(r,:));
    
    figure(figAllRois),subplot(1,numel(roiNames),r)
    bar([1,2],[acc_FF_mean(r),acc_BB_mean(r)],'facecolor',[0.15 0.45 0.75]),hold on
    bar([3.5,4.5],[acc_FB_mean(r),acc_BF_mean(r)],'facecolor',[0.55 0.55 0.55])
    set(gca,'xtick',[1 2 3.5 4.5],'xticklabels',{'FF','BB','FB','BF'},'fontsize',14)
    title(roiNamesTrue{r},'fontsize',18)
    plot(0:0.01:5.5,0.5*ones(length([0:0.01:5.5]),1),'color',[0.5 0.5 0.5],'linestyle','--')
    ylim([0.4 0.6]),xlim([0 5.5])
    
    aaa=mean([acc_BF(r,:);acc_FB(r,:)]);
    [h,p,ci,stats] = ttest(aaa-0.5)
    
end, clear r k s
save('results_1000perm_50')
clear responsePatterns

