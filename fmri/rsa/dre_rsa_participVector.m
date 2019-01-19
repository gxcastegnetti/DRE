%% dre_rsa_participVector
% ~~~
% GX Castegnetti --- 2019

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_pw';
betaid       = 'rsa_pulse_ima';
clear RESTOREDEFAULTPATH_EXECUTED

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
addpath(genpath([dir.rsaCod,fs,'drtoolbox']))
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

% directory with betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

%% subjects
subs = [4 5 7:9 13:17 19 21 23 25:26 29:32 34 35 37 39 40:43 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,3)];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_3L(dir,subs,taskOrd,bData);

%% which mask?

roiNames = {'rsaVal_LG_10mm','rsaVal_ACC_10mm','rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm'};
roiNames = {'rsaCon_vmPFC_10mm'};

%% prewhiten activity in the mask
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        % subject name
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        % SPM file from 1st level analysis
        subjSPMFile = [dir.beta,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        % load indices for normal ordering
        fooDir = [dir.dre,filesep,'out',filesep,'fmri',filesep,'rsa'];
        load([fooDir,filesep,'toNormalOrder',filesep,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
        
        % whiten betas
        try
            B = noiseNormaliseBeta_roi(SPM,subjMaskFile);
        catch
            keyboard
        end
        
        % take only those corresponding to conditions
        if size(B,1) == 124
            B = real(B([1:24,31:54,61:84,91:114],:));
        elseif size(B,1) == 268
            B = real(B([1:60,67:126,133:192,199:258],:));
            B = B(toNormalOrder,:);
        end
        
        B_struct.(roiNames{r}).(subjName) = B';
        
        % construct RDM
        rdm = squareform(pdist(B,'correlation'));
        RDM_brain(r,s).RDM = rdm;
    end
end, clear rdm fooDir subjSPMFile s SPM betaid

%% FIRE vs BOAT
for s = 1:length(subs)
    
    disp(['sub#',num2str(subs(s))])
    
    % load indices for normal ordering
    fooDir = [dir.dre,filesep,'out',filesep,'fmri',filesep,'rsa'];
    load([fooDir,filesep,'toNormalOrder',filesep,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
    
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
    objVal = [objVal_F; objVal_B];
    objCon = [objCon_F; objCon_B];
    
    %         X_B = B(1:120,:);
    %         X_F = B(121:end,:);
    %         mdl_F = fitrlinear(X_F(~isnan(objVal_F),:),objVal_F(~isnan(objVal_F)));
    %         mdl_B = fitrlinear(X_B(~isnan(objVal_B),:),objVal_B(~isnan(objVal_B)));
    %         betas_F = mdl_F.Beta;
    %         betas_B = mdl_B.Beta;
    %         r_betas(s) = corr(betas_F,betas_B);
    
    
    mdl_val = stepwiselm(B(~isnan(objVal),:),objVal(~isnan(objVal)));
    mdl_con = stepwiselm(B(~isnan(objCon),:),objCon(~isnan(objCon)));
    betas_val = mdl_val.Beta;
    betas_con = mdl_con.Beta;
    r_betas(s) = corr(betas_val,betas_con,'type','kendall');
    
    % permutation test
%     nPerm = 1000;
%     for i = 1:nPerm
%         fooRand_val = randperm(numel(objVal));
%         fooRand_con = randperm(numel(objVal));
%         objVal_perm = objVal(fooRand_val);
%         objCon_perm = objVal(fooRand_con);
%         
%         mdl_val_perm = fitrlinear(B(~isnan(objVal_perm),:),objVal_perm(~isnan(objVal_perm)));
%         mdl_con_perm = fitrlinear(B(~isnan(objCon_perm),:),objCon_perm(~isnan(objCon_perm)));
%         
%         betas_val_perm = mdl_val_perm.Beta;
%         betas_con_perm = mdl_con_perm.Beta;
%         r_betas_perm(s,i) = corr(betas_val_perm,betas_con_perm,'type','spearman');
%     end
end

% stats with real indices
[h,p,~,stats] = ttest(r_betas);
tVal_real = stats.tstat;

% for i = 1:nPerm
%     [h,p,~,stats] = ttest(r_betas_perm(:,i));
%     tVal_perm(i) = stats.tstat;
% end
% tVal_larger = sum(tVal_perm >= tVal_real);
% pVal = tVal_larger/nPerm;
