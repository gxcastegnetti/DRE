%% dre_classValue
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'class_v0';

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
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);

%% load response patterns and apply mask
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

% roiNames = {'box_w-16_16_16-0_-60_26'};
roiNames = {'HPC'};

% apply two masks: one for grey matter, one for ROI
for r = 1:length(roiNames)
    for s = 1:length(subs)
        subjName = ['SF',num2str(subs(s),'%03d')];
        roiMaskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        roiMask = spm_read_vols(spm_vol(roiMaskFile));
        
        % vectorise it
        roiMask = reshape(roiMask, 1, []);
        
        respPatt_foo = responsePatterns.(subjName)(logical(roiMask),:);
        respPatt.(['roi',num2str(r)]).(subjName) = respPatt_foo(~isnan(respPatt_foo(:,1)),:);
    end
end

subNames = fieldnames(respPatt.roi1);

for s = 1:length(subs)
    
    disp(['sub#',num2str(subs(s))])
    
    objIdx_1 = [bData(subs(s)).imagination(1).objIdx'; bData(subs(s)).imagination(3).objIdx'];
    objIdx_2 = [bData(subs(s)).imagination(2).objIdx'; bData(subs(s)).imagination(4).objIdx'];
    
    [~, objIdx_sort_1] = sort(objIdx_1);
    [~, objIdx_sort_2] = sort(objIdx_2);
    
    val_1 = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(3).val];
    val_2 = [bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(4).val];
    
    val_1 = val_1(objIdx_sort_1);
    val_2 = val_2(objIdx_sort_2);
    
    X = respPatt.roi1.(subNames{s})(:,1:120)';
    
    val_1(isnan(val_1)) = 25;
    val_2(isnan(val_2)) = 25;
    if taskOrd(s) == 1
        Y = val_1 + (0.00001*(1:120))';
    else
        Y = val_2 + (0.00001*(1:120))';
    end
    
    % separate into four levels
    Y_log = Y > median(Y);
    
    nTrials = 120;
    for lo = 1:nTrials
        
        testTrl = lo;
        
        X_test = X(testTrl,:);
        Y_test = Y_log(testTrl);
        X_train = X;
        X_train(testTrl,:) = [];
        Y_train = Y_log;
        Y_train(testTrl,:) = [];
        
        % B = fitlm(X,Y);
        [coeff, fitInfo] = lassoglm(X_train,double(Y_train),'binomial','alpha',1,'lambda',0.01);
        
        % Cau predictions
        foo = sum(coeff.*X_test(:)) + fitInfo.Intercept;
        prediction(testTrl) = round(1/(1+exp(-foo)));
        
    end
    
    performance(s) = sum(prediction(:) == Y_log)/nTrials;
    
end
