%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_pw';
betaid       = 'rsa_pulse_ons0';

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
subs = [5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,9),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% define betas and apply mask
roiNames = {'l_hpc'};
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
        B = real(B([1:60,67:126,133:192,199:258],:));
        B = B(toNormalOrder,:);
%         B = B(ordData(subs(s)).norm2val_cont,:);
        
        % construct RDM
        rdm = squareform(pdist(B,'euclidean'));
        RDM_struct(s).RDM = rdm;
        RDM_struct(s).name = ['sub#',num2str(subs(s),'%03d')];
        RDM_struct(s).color = [0 0 1];
        RDM_all(:,:,s) = rdm;
    end
end
% figureRDMs(RDM_struct,userOptions)

RDM_mean.RDM = mean(RDM_all,3);
RDM_mean.name = 'Mahalanobis';
RDM_mean.color = [0 0 1];
figureRDMs(RDM_mean,userOptions)

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