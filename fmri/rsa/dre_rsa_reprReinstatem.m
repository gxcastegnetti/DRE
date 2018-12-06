%% dre_reprReinstatem
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'reprReinstatement_v0';

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

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);

%% ROI
roiNames = {'lingual','l_hpc','r_hpc','pcc','mcc','pfc_vm','ofc'};
% roiNames = {'lingual','itc','phpc','l_hpc','r_hpc','angular','par_inf','ins_la','pcc','mcc','acc','caudate','putamen','subgenual','pfc_vm','ofc'};
roiNames = {'lp_hpc','la_hpc','rp_hpc','ra_hpc'};

%% betas for imagination and choice
dir.betaIma = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,'rsa_pulse_ima',fs,'none'];
dir.betaCho = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,'rsa_pulse_choice',fs,'none'];

%% apply two masks: one for grey matter, one for ROI
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        % subject name
        subjName = ['SF',num2str(subs(s),'%03d')];
        
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
        
        %%%%%%%%%%%%%%%
        % imagination %
        %%%%%%%%%%%%%%%
        
        % SPM file from 1st level analysis
        subjSPMFile_ima = [dir.betaIma,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile_ima)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        % whiten betas
        B_ima_foo = noiseNormaliseBeta_roi(SPM,subjMaskFile);
        
        % sort them
        B_ima_foo = real(B_ima_foo([1:60,67:126,133:192,199:258],:));
        B_ima_foo = B_ima_foo(toNormalOrder,:);
        B_ima.(roiNames{r}).(subjName) = B_ima_foo';
        clear subjMaskFile SPM B_ima_foo
        
        %%%%%%%%%%
        % choice %
        %%%%%%%%%%
        
        % SPM file from 1st level analysis
        subjSPMFile_cho = [dir.betaCho,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile_cho)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        % whiten betas
        B_cho_foo = noiseNormaliseBeta_roi(SPM,subjMaskFile);
        
        % sort them
        B_cho_foo = real(B_cho_foo([1:24,31:54,61:84,91:114],:));
        B_cho.(roiNames{r}).(subjName) = B_cho_foo';
        clear subjMaskFile SPM B_cho_foo
    end
end

%% define field names
roiNames = fieldnames(B_cho);
subNames = fieldnames(B_cho.(roiNames{1}));

%% correlations between response pattern during choice and corresponding imaginations
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        if taskOrd(s) == 1
            trlOffset = [zeros(24,1);ones(24,1);zeros(24,1);ones(24,1)];
        else
            trlOffset = [ones(24,1);zeros(24,1);ones(24,1);zeros(24,1)];
        end
        
        % stack objects presented in choices during all sessions
        ObjPresentedAll = [bData(subs(s)).choice(1).objPresented;
            bData(subs(s)).choice(2).objPresented;
            bData(subs(s)).choice(3).objPresented;
            bData(subs(s)).choice(4).objPresented];
        
        % stack which object was chosen
        ObjChosenAll = [bData(subs(s)).choice(1).objChosen;
            bData(subs(s)).choice(2).objChosen;
            bData(subs(s)).choice(3).objChosen;
            bData(subs(s)).choice(4).objChosen];
        
        for trlCho = 1:96
            % take response pattern during choice
            respPattCho_trl = B_cho.(roiNames{r}).(subNames{s})(:,trlCho);
            
            % take corresponding response patterns during imagination
            obj_1_idx = ObjPresentedAll(trlCho,1) + 120*trlOffset(trlCho);
            obj_2_idx = ObjPresentedAll(trlCho,2) + 120*trlOffset(trlCho);
            respPattIma_trl_1 = B_ima.(roiNames{r}).(subNames{s})(:,obj_1_idx);
            respPattIma_trl_2 = B_ima.(roiNames{r}).(subNames{s})(:,obj_2_idx);
            if ObjChosenAll(trlCho) == -1
                respPattIma_trl_chosen = respPattIma_trl_1;
                respPattIma_trl_unchos = respPattIma_trl_2;
            elseif ObjChosenAll(trlCho) == 1
                respPattIma_trl_chosen = respPattIma_trl_2;
                respPattIma_trl_unchos = respPattIma_trl_1;
            end
            
            % compute correlation with chosen and unchosen
            corr_choiceVSchosen(trlCho) = corr(respPattCho_trl,respPattIma_trl_chosen);
            corr_choiceVSunchos(trlCho) = corr(respPattCho_trl,respPattIma_trl_unchos);
        end
        %         corr_choiceVSchosen_sub(r,s) = mean(corr_choiceVSchosen);
        %         corr_choiceVSunchos_sub(r,s) = mean(corr_choiceVSunchos);
        corr_choiceVdiffere_sub(r,s) = mean(corr_choiceVSchosen - corr_choiceVSunchos);
    end
end

%% ttest
for r = 1:length(roiNames)
    scores = corr_choiceVdiffere_sub(r,:);
    [h,p(r),~,~] = ttest(scores,0,'Tail','right');
end
