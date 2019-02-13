%% dre_reprReinstatemVsConf
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
subs = [4 5 7 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39:43 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,4)];

subsBest = sort([24,19,6,4,22,12,11,21,18,30,25,1,16,23,2]);
subsWors = sort([15,7,5,8,3,10,28,13,17,20,14,9,26,29,27]);

% subs = subs(subsBest);
% taskOrd = taskOrd(subsBest);

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);

%% ROI
roiNames = {'lp_hpc','la_hpc','rp_hpc','ra_hpc'};
%
% roiNames = {'l_hpc','r_hpc'};
% roiNames = {'rsaOid_Occ_10mm','rsaVal_ITG'};
% roiNames = {'l_phpc','r_phpc'};


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
        ObjPresented_allSession = [bData(subs(s)).choice(1).objPresented;
            bData(subs(s)).choice(2).objPresented;
            bData(subs(s)).choice(3).objPresented;
            bData(subs(s)).choice(4).objPresented];
        
        ObjPresented_allSession = ObjPresented_allSession + 120*trlOffset;
        
        % stack which object was chosen
        ObjChosen_allSessions = [bData(subs(s)).choice(1).objChosen;
            bData(subs(s)).choice(2).objChosen;
            bData(subs(s)).choice(3).objChosen;
            bData(subs(s)).choice(4).objChosen];
        
        corr_choiceVSleft = nan(96,1);
        corr_choiceVSrigh = nan(96,1);
        for trlCho = 1:96
            % take response pattern during choice
            respPattCho_trl = B_cho.(roiNames{r}).(subNames{s})(:,trlCho);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % take corresponding response patterns during imagination %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % take indices corresponding to the two objects proposed
            obj_1_idx = ObjPresented_allSession(trlCho,1);
            obj_2_idx = ObjPresented_allSession(trlCho,2);
            
            % take the corresponding activation pattern during imagination
            respPattIma_leftSide = B_ima.(roiNames{r}).(subNames{s})(:,obj_1_idx);
            respPattIma_righSide = B_ima.(roiNames{r}).(subNames{s})(:,obj_2_idx);

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % uncomment to compare H/L value instead of chosen/unchosen %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % load indices for normal ordering
            dirdre = '/Users/gcastegnetti/Desktop/stds/DRE';
            fooDir = [dirdre,filesep,'out',filesep,'fmri',filesep,'rsa'];
            load([fooDir,filesep,'toNormalOrder',filesep,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
            
            % stack which object was chosen
            val_allSessions = [bData(subs(s)).imagination(1).con;
                bData(subs(s)).imagination(2).con;
                bData(subs(s)).imagination(3).con;
                bData(subs(s)).imagination(4).con];
            
            val_allSessions = val_allSessions(toNormalOrder);
            
            % take value of the two objects
            val_1 = val_allSessions(obj_1_idx);
            val_2 = val_allSessions(obj_2_idx);
            
            valDiff(trlCho) = abs(val_1);
            
            % compute correlation with chosen and unchosen
            corr_choiceVSleft(trlCho) = corr(respPattCho_trl,respPattIma_leftSide);
            corr_choiceVSrigh(trlCho) = corr(respPattCho_trl,respPattIma_righSide);
        end
        %         corr_choiceVSchosen_sub(r,s) = mean(corr_choiceVSchosen);
        %         corr_choiceVSunchos_sub(r,s) = mean(corr_choiceVSunchos);
        corr_choiceVSboth = nanmean([corr_choiceVSleft],2);
        corr_reinstatVsConf(r,s) = corr(corr_choiceVSboth(:),valDiff(:),'type','spearman','rows','complete');
    end
end

%% chosen - unchosen

% t-test
for r = 1:length(roiNames)
    scores = corr_reinstatVsConf(r,:);
    %     [h,p(r),~,~] = ttest(scores,0,'Tail','right');
    [p(r),h,stats] = signrank(scores,0,'Tail','right');
end

% plot
figure('color',[1 1 1])
% roiNamesTrue = {'MidOcc','Lingual','Lingual (uni)','lp HPC','la HPC','rp HPC','ra HPC','PCC','MCC','ACC','OFC'};
roiNamesTrue = {'LP','LA','RP','RA'};
corrMean = mean(corr_choiceVSboth,2);
corrSem = std(corr_choiceVSboth,0,2)/sqrt(numel(subs));
bar(1:numel(roiNames),mean(corr_choiceVSboth,2),0.4), hold on
errorbar(1:numel(roiNames),corrMean,corrSem,'k.')
set(gca,'fontsize',18,'xtick',1:numel(roiNames),'xticklabels',roiNamesTrue)
ylabel('r(choice,chosObj)-r(choice,unchObj)')


