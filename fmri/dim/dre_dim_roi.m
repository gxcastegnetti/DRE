%% dre_dim_roi
% ~~~
% performs roi dimensionality analysis
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'dim_roi_test_1';
rsaAnalysisName = 'rsa_roi_pulse_ons0';

%% folders
fs      = filesep;
dir.dim = pwd;
idcs    = strfind(dir.dim,'/');
dir.dre = dir.dim(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsa = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'dim',fs,'roi'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath([dir.sta,fs,'routines']))
addpath(genpath('functionalDimensionality'))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% extract behaviour
bData = dre_extractData(dir,subs,taskOrd,0);

%% load response patterns computed in dre_rsa_roi_run
load([dir.rsaOut,fs,rsaAnalysisName,fs,'rsaPatterns_roi.mat'])
respPatt = responsePatterns; clear responsePatterns

%% find roi and sub names
roiNames = fieldnames(respPatt);
subNames = fieldnames(respPatt.(roiNames{1}));

%% rearrange response patterns in presentation order
for sub = 1:length(subNames)
    for roi = 1:length(roiNames)
        
        % reorder objects according to presentation
        sessType1 = bData(subs(sub)).sessType{1};
        idx_s1 = sort(bData(subs(sub)).imagination(1).objIdx);
        if strcmp(sessType1,'B')
            idx_s1 = idx_s1 + 120;
        end
        
        sessType2 = bData(subs(sub)).sessType{2};
        idx_s2 = sort(bData(subs(sub)).imagination(2).objIdx);
        if strcmp(sessType2,'B')
            idx_s2 = idx_s2 + 120;
        end
        
        sessType3 = bData(subs(sub)).sessType{3};
        idx_s3 = sort(bData(subs(sub)).imagination(3).objIdx);
        if strcmp(sessType3,'B')
            idx_s3 = idx_s3 + 120;
        end
        
        sessType4 = bData(subs(sub)).sessType{4};
        idx_s4 = sort(bData(subs(sub)).imagination(4).objIdx);
        if strcmp(sessType4,'B')
            idx_s4 = idx_s4 + 120;
        end
        
        % for every subject, this should be nVox x nCond x nRuns
        respPatt_presentOrder{roi}{sub}{1} = respPatt.(roiNames{roi}).(subNames{sub})(:,idx_s1); %#ok<*SAGROW>
        respPatt_presentOrder{roi}{sub}{2} = respPatt.(roiNames{roi}).(subNames{sub})(:,idx_s2);
        respPatt_presentOrder{roi}{sub}{3} = respPatt.(roiNames{roi}).(subNames{sub})(:,idx_s3);
        respPatt_presentOrder{roi}{sub}{4} = respPatt.(roiNames{roi}).(subNames{sub})(:,idx_s4);
    end
end

%% PCA
explVar_threshold = 75;
for roi = 1:length(roiNames)
    for sub = 1:length(subNames)
        for sess = 1:4
            [~,~,~,~,explVar_vec,~] = pca(respPatt_presentOrder{roi}{sub}{sess});
            for i = 1:length(explVar_vec)
                explVar_soFar = sum(explVar_vec(1:i));
                if explVar_soFar > explVar_threshold
                    break
                end
            end
            explVar_sess{roi}(sub,sess) = i;
        end
    end
    explVar_sub{roi} = mean(explVar_sess{roi},2);
end

% t-test
tTestResults = zeros(length(roiNames));
for roi_1 = 1:length(roiNames)
    for roi_2 = 1:length(roiNames)
        [h,p,ci,stats] = ttest(explVar_sub{roi_1},explVar_sub{roi_2});
        if h == 1
            tTestResults(roi_1,roi_2) = stats.tstat;
        end
    end
end

keyboard

%
% % load gm mask
% binaryMask = [dir.msk,fs,'rgm.nii'];
% for r = 1:length(roiNames)
%     [bestn{r},r_outer{r}, r_alter{r}, test_tfce{r}] = functional_dimensionality(respPatt_presentOrder{r}, 0);
% end
