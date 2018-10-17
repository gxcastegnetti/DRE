%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_sphere';
dirBeta = 'rsa_pulse_ons0';

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
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% mask brains
roiNames = {'left_post_hpc_9'};

filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')


for r = 1:length(roiNames)
    for s = 1:length(subs)
        subjName = ['SF',num2str(subs(s),'%03d')];
        maskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        mask = spm_read_vols(spm_vol(maskFile));
        
        % vectorise it
        mask = reshape(mask, 1, []);
        mask = logical(mask');
        
        respPatt{s} = responsePatterns.(subjName)(mask,:);
    end
end


%% construct RDMs
% RDMs_data = constructRDMs(respPatt, 'SPM', userOptions);

for s = 1:28
    RDMs_data{s} = squareform(pdist(respPatt{s}', 'correlation'));
end

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];
for s = 1:length(subs)
    RDMs_model(1,s).name = 'val';
    RDMs_model(1,s).RDM = RDMs_models{s}.val; %#ok<*SAGROW>
    RDMs_model(1,s).color = [0 1 0];
    RDMs_model(2,s).name = 'fam';
    RDMs_model(2,s).RDM = RDMs_models{s}.fam;
    RDMs_model(2,s).color = [0 1 0];
    RDMs_model(3,s).name = 'oid';
    RDMs_model(3,s).RDM = 1-mat_ID;
    RDMs_model(3,s).color = [0 1 0];
    RDMs_model(4,s).name = 'cxt';
    RDMs_model(4,s).RDM = RDMs_models{s}.cxt;
    RDMs_model(4,s).color = [0 1 0];
    
%     for m = 1:4
%         modelRDMs_ltv(m,:) = permute(unwrapRDMs(vectorizeRDMs(RDMs_model(m,s))), [3 2 1]);
%     end
end

%% for every region and sub, correlate RDM and model
scoreNames = {'val','fam','ID','cxt'};

for s = 1:size(RDMs_data,2)
    for m = 1:size(RDMs_model,1)
        a = vectorizeRDM(RDMs_data{s});
        b = permute(unwrapRDMs(vectorizeRDM(RDMs_model(m,s).RDM)),[3 2 1]);
        rL2(m,s) = corr(b',a','rows','pairwise','type','Spearman');
        rL2(m,s) = fisherTransform(rL2(m,s));
    end
end
figure
bar(mean(rL2(:,:,r),2)),set(gca,'xticklabel',scoreNames)
set(gca,'fontsize',16)
ylim([-0.01 0.01])

%% ttest
for m = 1:size(RDMs_model,1)
    
    scores = rL2(m,:,r);
    [h,p(m),~,~] = ttest(scores);
end

%%





