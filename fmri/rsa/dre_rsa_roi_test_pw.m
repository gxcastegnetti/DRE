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
subs = [4 5 7:9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% define betas and apply mask
roiNames = {'box_w-16_16_16-0_-60_26'};
roiNames = {'thalamus'};
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

% apply two masks: one for grey matter, one for ROI
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
    subjMaskFile.fname = [dir.mskOut,fs,roiNames{1},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{1},'.nii'];
    
    % whiten betas
    B = noiseNormaliseBeta_roi(SPM,subjMaskFile);
    
    % take only those corresponding to conditions
    B = real(B([1:60,67:126,133:192,199:258],:));
    B = B(toNormalOrder,:);
    
    fooDir = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
    save([fooDir,fs,'toNormalOrder',fs,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
    %         B = B(ordData(subs(s)).norm2val_cont,:);
    
    % construct RDM
    rdm = squareform(pdist(B,'correlation'));
%     RDM_struct(s).RDM = rdm;
%     RDM_struct(s).name = ['sub#',num2str(subs(s),'%03d')];
%     RDM_struct(s).color = [0 0 1];
    RDM_brain(s).RDM = rdm;
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

%% for every region and sub, correlate RDM and model
scoreNames = {'val','con','fam','ID','cxt'};
for s = 1:size(RDM_brain,2)
    for m = 1:size(RDM_model,1)
        a = vectorizeRDM(RDM_brain(s).RDM);
        b = vectorizeRDM(RDM_model(m,s).RDM);
        corrRoiModel(m,s) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
        %         rL2(m,s) = fisherTransform(rL2(m,s));
    end
end

%% ttest
for m = 1:size(RDM_model,1)
    scores = corrRoiModel(m,:);
    [h,p(m),~,~] = ttest(scores);
end