%% dre_rsa_roi_test_hpc_change
% ~~~
% GX Castegnetti --- 2018

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
addpath([dir.rsaCod,fs,'routines']),
addpath(genpath([dir.rsaCod,fs,'drtoolbox']))
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

% directory with betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

%% subjects
subs = [4 5 7:9 13:17 19 21 23 25:26 29:32 34 35 37 39 40:43 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,4)];

subsBest = sort([24,19,6,4,22,12,11,21,18,30,25,1,16,23,2]);
subsWors = sort([15,7,5,8,3,10,28,13,17,20,14,9,26,29,27]);

slopes = [0.244644048172450,0.223946840860126,0.208531637926186,0.400781101520874,0.214651122900286,...
    0.527213302443985,0.222262713073930,0.214265149192790,0.113871322525068,0.208145500659497,...
    0.303896467054236,0.310876174215477,0.157495707488989,0.123082743374363,0.223259451529286,...
    0.233042188472160,0.134694830218049,0.278532556059753,0.548450276749430,0.128458679610595,...
    0.286603706120169,0.333069644635516,0.227243248191941,0.650551192022550,0.262801817706755,...
    0.150765971850011,0.0835213995036290,0.163443588603048,0.192293815632133,0.276329355358742];

% subs = subs(subsBest);
% taskOrd = taskOrd(subsBest);

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);

%% which mask?
roiNames = {'rsaVal_ACC_10mm','rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm'};
roiNames = {'lp_hpc','rp_hpc'};

%% prewhiten activity in the mask
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
        toNormalOrder(s,:) = [objIdx_F,objIdx_B];
        
        %         fooDir = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
        %         save([fooDir,fs,'toNormalOrder',fs,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
        
        % SPM file from 1st level analysis
        subjSPMFile = [dir.beta,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        %         % subjective gm mask
        %         subjGreyFile.fname = [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
        %
        %         % take intersection between the two
        
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
            B = B(toNormalOrder(s,:),:);
        end
        
        B_struct.(roiNames{r}).(subjName) = B';
        
        % construct RDM
        rdm = squareform(pdist(B,'correlation'));
        RDM_brain(r,s).RDM = rdm;
    end
end, clear rdm fooDir subjSPMFile objIdx_F objIdx_B obj objIdx_sess_F objIdx_sess_B s SPM betaid

%% ROI-model correlations

%%%%%%%%%%%%%%%%%
% create models %
%%%%%%%%%%%%%%%%%
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);

% % for every region and sub, correlate RDM and model
% for r = 1:length(roiNames)
%     for s = 1:size(RDM_brain,2)
%         
%         % take only similarities between identical objects
%         brain_offDiagMat = RDM_brain(r,s).RDM(1:120,121:end);
%         brain_dissSameItem = diag(brain_offDiagMat);
%         
%         model_offDiagMat = RDMs_models{s}.val(1:120,121:end);
%         model_dissSameItem = diag(model_offDiagMat);
%         
%         corrValChange(r,s) = corr(brain_dissSameItem(:),model_dissSameItem(:),'rows','pairwise','type','spearman');
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ttest of correlations for each ROI and model %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for r = 1:length(roiNames)
%     
%     scores = corrValChange(r,:);
%     [pCorr(r),h,~] = signrank(scores,0,'Tail','left');
%     [rPerf.r(r),rPerf.p(r)] = corr(scores(:),slopes(:),'type','spearman');
%     
% end, clear r m mat_ID

% test orthogonalisation in the hippocampus
for r = 1:length(roiNames)
    for s = 1:numel(subs)
        
        % take only similarities between identical objects
        geom_F = vectorizeRDM(RDM_brain(r,s).RDM(1:120,1:120));
        geom_B = vectorizeRDM(RDM_brain(r,s).RDM(121:end,121:end));       
        
        corrValChange(r,s) = corr(geom_F(:),geom_B(:),'rows','pairwise','type','spearman');
    end
end

for r = 1:length(roiNames)
    
    scores = corrValChange(r,:);
    [pCorr(r),h,~] = signrank(scores,0,'Tail','left');
    [rPerf.r(r),rPerf.p(r)] = corr(scores(:),slopes(:),'type','spearman');
    
end, clear r m mat_ID

