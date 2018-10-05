%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_test_1';

%% folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.msk  = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];
addpath([dir.root,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% load response patters computed in dre_rsa_roi_run
load([dir.out,fs,analysisName,fs,'rsaPatterns_roi.mat'],'responsePatterns')

%% construct RDMs
RDMs_data = constructRDMs(responsePatterns, 'SPM', userOptions);
RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');

%% plot RDMs
% matrices
figureRDMs(RDM_average,userOptions)

% dendrograms
% dendrogramConditions(RDM_average,userOptions)

% MDS
% MDSConditions(RDM_average,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];
mat_cxt = [zeros(120), ones(120);
    ones(120), zeros(120)];
for s = 1:length(subs)
    RDMs_model(1,s).name = 'value';
    RDMs_model(1,s).RDM = RDMs_models{s}.val; %#ok<*SAGROW>
    RDMs_model(1,s).color = [0 1 0];
    RDMs_model(2,s).name = 'confidence';
    RDMs_model(2,s).RDM = RDMs_models{s}.con;
    RDMs_model(2,s).color = [0 1 0];
    RDMs_model(3,s).name = 'familiarity';
    RDMs_model(3,s).RDM = RDMs_models{s}.fam;
    RDMs_model(3,s).color = [0 1 0];
    RDMs_model(4,s).name = 'price';
    RDMs_model(4,s).RDM = RDMs_models{s}.pri;
    RDMs_model(4,s).color = [0 1 0];
    RDMs_model(5,s).name = 'obj ID';
    RDMs_model(5,s).RDM = 1-mat_ID;
    RDMs_model(5,s).color = [0 1 0];
    RDMs_model(6,s).name = 'context';
    RDMs_model(6,s).RDM = mat_cxt;
    RDMs_model(6,s).color = [0 1 0];
end


%% for every region and sub, correlate RDM and model
regionNames = {'hpc','mpfc','lingual'};
scoreNames = {'val','con','fam','pri','ID','cxt'};
h{1} = figure('color',[1 1 1]);
h{2} = figure('color',[1 1 1]);
h{3} = figure('color',[1 1 1]);
for r = 1:size(RDMs_data,1)
    for s = 1:size(RDMs_data,2)
        for m = 1:size(RDMs_model,1)
            a = vectorizeRDM(RDMs_data(r,s).RDM);
            b = vectorizeRDM(RDMs_model(m,s).RDM);
            rL2(m,s) = corr(a',b','rows','complete','type','Spearman');
        end
        figure(h{r})
        subplot(5,6,s),bar(rL2(:,s)),set(gca,'xticklabel',scoreNames)
        title([regionNames{r},' - sub#',num2str(subs(s),'%03d')])
    end
    figure('color',[1 1 1]),bar(mean(rL2,2)),set(gca,'xticklabel',scoreNames),title(regionNames{r})
end






keyboard
for m = 1:size(RDMs_data,1)
    for s = 1:length(subs)
        
        RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_con{s});
        corrMat = RDMCorrMat(RDMs, 1);
        c(m,s) = corrMat(1,2);
        
        %         RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_val{s});
        %         corrMat = RDMCorrMat(RDMs, 1);
        %         c(m,s) = corrMat(1,2);
        %
        %         RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_val{s});
        %         corrMat = RDMCorrMat(RDMs, 1);
        %         c(m,s) = corrMat(1,2);
        %
    end
    [h,p(m),~,~] = ttest(c(m,:))
end

%% test goal
% find ROIs considered in this analysis
roiNames = fieldnames(responsePatterns);

for i = 1:length(roiNames)
    % for each subject take the vector of presented items
    for s = 1:length(subs)
        
        % find indices of items presented in each session. Add 120 to the boat
        % conditions because in responsePatters the 120 boat activity
        % patterns are stacked after the 120 fire patterns.
        if taskOrd(s) == 1
            idx_S1 = bData(subs(s)).imagination(1).fire.objIdx;
            idx_S2 = bData(subs(s)).imagination(2).boat.objIdx + 120;
            idx_S3 = bData(subs(s)).imagination(3).fire.objIdx;
            idx_S4 = bData(subs(s)).imagination(4).boat.objIdx + 120;
        else
            idx_S1 = bData(subs(s)).imagination(1).boat.objIdx + 120;
            idx_S2 = bData(subs(s)).imagination(2).fire.objIdx;
            idx_S3 = bData(subs(s)).imagination(3).boat.objIdx + 120;
            idx_S4 = bData(subs(s)).imagination(4).fire.objIdx;
        end
        
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        % fetch activation patterns
        actPatt_S1 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S1);
        actPatt_S2 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S2);
        actPatt_S3 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S3);
        actPatt_S4 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S4);
        
        corr_S1S2 = corr(actPatt_S1,actPatt_S2);
        corr_S1S3 = corr(actPatt_S1,actPatt_S3);
        foo(s) = mean(mean(corr_S1S3)) - mean(mean(corr_S1S2));
        
    end
    
    roiNames{i}
    [h,p,~,~] = ttest(foo)
    
end

