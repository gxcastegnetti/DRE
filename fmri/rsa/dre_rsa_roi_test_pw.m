%% dre_rsa_roi_test_pw
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
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'drtoolbox']))
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

% directory with betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

%% subjects
subs = [4 5 7:9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_3L(dir,subs,taskOrd,bData);

%% which mask?
% roiNames = {'l_ling','lp_itc','lp_hpc','rp_hpc','la_hpc','ra_hpc','mcc','sma','rp_ins','la_ins','l_dlpfc','r_dlpfc','l_ofc','vmpfc_ima'};
% roiNames = {'l_ling','lp_itc','l_hpc','r_hpc','lp_hpc','rp_hpc','mcc','sma','rp_ins','la_ins','l_dlpfc','r_dlpfc','vmpfc_ima_p','l_ofc','ofc_conf'};
roiNames = {'l_ling','lp_hpc','rp_hpc','mcc','rp_ins','vmpfc_ima_p','ofc_conf','l_ofc'};

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
        
        % SPM file from 1st level analysis
        subjSPMFile = [dir.beta,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
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
        
        % fooDir = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
        % save([fooDir,fs,'toNormalOrder',fs,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
        
        % construct RDM
        rdm = squareform(pdist(B,'correlation'));
        RDM_brain(r,s).RDM = rdm;
    end
end, clear rdm fooDir subjSPMFile objIdx_F objIdx_B obj objIdx_sess_F objIdx_sess_B s SPM betaid

%% plot ROI RDM
% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

% average RDMs across subjects and plot
figure('color',[1,1,1])
for r = 1:numel(roiNames)
    
    % compute ROI-specific average
    foo = zeros(size(B,1),size(B,1));
    for s = 1:numel(subs)
        newOrd = ordData(subs(s)).norm2val_cont;
        foo = foo + RDM_brain(r,s).RDM(newOrd,newOrd)/numel(subs);
    end
    
    % fill struct for the RSA toolbox
    RDM_brain_avg(r).color = [0 0 1];
    RDM_brain_avg(r).name = roiNames{r};
    RDM_brain_avg(r).RDM = foo;
    
    % plot
    subplot(3,5,r)
    RDM_brain_avg(r).RDM(RDM_brain_avg(r).RDM==0)=nan;
    foo = nanmean(RDM_brain_avg(r).RDM(:));
    imagesc(RDM_brain_avg(r).RDM),title(roiNames{r}),caxis([foo-1 foo+1])
end


%% ROI-model correlations

%%%%%%%%%%%%%%%%%
% create models %
%%%%%%%%%%%%%%%%%
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];

%%%%%%%%%%%%%%%%%%
% if imagination %
%%%%%%%%%%%%%%%%%%
if size(RDM_brain(1,1).RDM,1) == 240
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
    
    % for every region and sub, correlate RDM and model
    scoreNames = {'Value','Confidence','Familiarity','Object ID','Context'};
    for r = 1:length(roiNames)
        for s = 1:size(RDM_brain,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
                % rL2(m,s) = fisherTransform(rL2(m,s));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%
    % else, if choice %
    %%%%%%%%%%%%%%%%%%%
elseif size(RDM_brain(1,1).RDM,1) == 96
    for s = 1:length(subs)
        RDM_model(1,s).name = 'dVal';
        RDM_model(1,s).RDM = RDMs_models{s}.choice.dVal; %#ok<*SAGROW>
        RDM_model(1,s).color = [0 1 0];
        RDM_model(2,s).name = 'cMun';
        RDM_model(2,s).RDM = RDMs_models{s}.choice.cMun;
        RDM_model(2,s).color = [0 1 0];
        RDM_model(3,s).name = 'chos';
        RDM_model(3,s).RDM = RDMs_models{s}.choice.Chos;
        RDM_model(3,s).color = [0 1 0];
        RDM_model(4,s).name = 'unch';
        RDM_model(4,s).RDM = RDMs_models{s}.choice.Unch;
        RDM_model(4,s).color = [0 1 0];
        RDM_model(5,s).name = 'ccxt';
        RDM_model(5,s).RDM = RDMs_models{s}.choice.ccxt;
        RDM_model(5,s).color = [0 1 0];
    end
    
    % for every region and sub, correlate RDM and model
    scoreNames = {'Value diff.','Chosen - unch.','Chosen value','Unch. value','Context'};
    for r = 1:length(roiNames)
        for s = 1:size(RDM_brain,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','Spearman');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ttest of correlations for each ROI and model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(roiNames)
    for m = 1:size(RDM_model,1)
        scores = corrRoiModel(r,:,m);
        [h,pCorr(r,m),~,~] = ttest(scores,0,'Tail','right');
    end
end, clear r m mat_ID

%%%%%%%%%%%%%%%%%%%%%
% plot mean and sem %
%%%%%%%%%%%%%%%%%%%%%
means = squeeze(mean(corrRoiModel,2));
sems  = squeeze(std(corrRoiModel,0,2)/sqrt(numel(subs)));

% means = means(end-1:end,1:2);
% sems = sems(end-1:end,1:2);

roiNamesTrue = roiNames;
% roiNamesTrue = roiNames;
figure('color',[1 1 1])
hb = bar(means); hold on
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData',means(:,ib),sems(:,ib),'k.')
    %     hb(ib).FaceColor = myColors(ib,:);
end, clear ib
legend({'Value','Confid.'},'location','northwest'),set(gca,'fontsize',18,'xtick',1:numel(roiNamesTrue),...
    'xticklabels',roiNames)
ylabel('Correlation(ROI, model)')


%% ROI-ROI comparison in terms of correlation with mdoel
for r1 = 1:length(roiNames)
    for r2 = 1:length(roiNames)
        corrRoiVal_r1 = squeeze(corrRoiModel(r1,:,1));
        corrRoiVal_r2 = squeeze(corrRoiModel(r2,:,1));
        [~,p(r1,r2),~,stats] = ttest(corrRoiVal_r2 - corrRoiVal_r1,0,'tail','right');
        tValues(r1,r2) = stats.tstat;
    end
end
isSignificant = p < 0.05;
tValues = tValues.*isSignificant;
figure('color',[1 1 1]),imagesc(tValues)


%% ROI-ROI correlations

whatScore = 'value';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take trials with hi/low val, con, fam %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(subs)
    
    % extract scores
    objVal = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(2).val; bData(subs(s)).imagination(3).val; bData(subs(s)).imagination(4).val];
    objCon = [bData(subs(s)).imagination(1).con; bData(subs(s)).imagination(2).con; bData(subs(s)).imagination(3).con; bData(subs(s)).imagination(4).con];
    objFam = [bData(subs(s)).imagination(1).fam; bData(subs(s)).imagination(2).fam; bData(subs(s)).imagination(3).fam; bData(subs(s)).imagination(4).fam];
    
    %     objVal = [bData(subs(s)).choice(1).valCho; bData(subs(s)).choice(2).valCho; bData(subs(s)).choice(3).valCho; bData(subs(s)).choice(4).valCho];
    
    % fix nans
    objVal(isnan(objVal)) = ceil(50*rand);
    objCon(isnan(objCon)) = ceil(50*rand);
    objFam(isnan(objFam)) = ceil(50*rand);
    
    % add a constant for univoque determination of percentile
    Y_val = objVal + (0.000001*(1:numel(objVal)))';
    Y_con = objCon + (0.000001*(1:numel(objCon)))';
    Y_fam = objFam + (0.000001*(1:numel(objFam)))';
    
    % sort them
    Y_val = Y_val(toNormalOrder(s,:));
    Y_con = Y_con(toNormalOrder(s,:));
    Y_fam = Y_fam(toNormalOrder(s,:));
    
    % find percentiles
    pctile_val_33 = prctile(Y_val,100/3);
    pctile_val_66 = prctile(Y_val,200/3);
    pctile_con_33 = prctile(Y_con,100/3);
    pctile_con_66 = prctile(Y_con,200/3);
    pctile_fam_33 = prctile(Y_fam,100/3);
    pctile_fam_66 = prctile(Y_fam,200/3);
    
    % separate trials according to value
    val_LO(s,:) = Y_val < pctile_val_33;
    val_HI(s,:) = Y_val > pctile_val_66;
    
    % separate trials according to value
    con_LO(s,:) = Y_con < pctile_con_33;
    con_HI(s,:) = Y_con > pctile_con_66;
    
    % separate trials according to value
    fam_LO(s,:) = Y_fam < pctile_fam_33;
    fam_HI(s,:) = Y_fam > pctile_fam_66;
    
end


for r1 = 1:length(roiNames)
    for r2 = 1:length(roiNames)
        for s = 1:size(RDM_brain,2)
            
            % subject name
            subjName = ['SF',num2str(subs(s),'%03d')];
            
            % take activity in the two ROIs
            if strcmp(whatScore,'value')
                roi1_HI = B_struct.(roiNames{r1}).(subjName)(:,val_HI(s,:));
                roi1_LO = B_struct.(roiNames{r1}).(subjName)(:,val_LO(s,:));
                roi2_HI = B_struct.(roiNames{r2}).(subjName)(:,val_HI(s,:));
                roi2_LO = B_struct.(roiNames{r2}).(subjName)(:,val_LO(s,:));
            elseif strcmp(whatScore,'confidence')
                roi1_HI = B_struct.(roiNames{r1}).(subjName)(:,con_HI(s,:));
                roi1_LO = B_struct.(roiNames{r1}).(subjName)(:,con_LO(s,:));
                roi2_HI = B_struct.(roiNames{r2}).(subjName)(:,con_HI(s,:));
                roi2_LO = B_struct.(roiNames{r2}).(subjName)(:,con_LO(s,:));
            elseif strcmp(whatScore,'familiarity')
                roi1_HI = B_struct.(roiNames{r1}).(subjName)(:,fam_HI(s,:));
                roi1_LO = B_struct.(roiNames{r1}).(subjName)(:,fam_LO(s,:));
                roi2_HI = B_struct.(roiNames{r2}).(subjName)(:,fam_HI(s,:));
                roi2_LO = B_struct.(roiNames{r2}).(subjName)(:,fam_LO(s,:));
            end
            
            % now create RDMs
            rdm_roi1_HI = pdist(roi1_HI','correlation');
            rdm_roi1_LO = pdist(roi1_LO','correlation');
            rdm_roi2_HI = pdist(roi2_HI','correlation');
            rdm_roi2_LO = pdist(roi2_LO','correlation');
%             corrRoiRoi_HI(r1,r2,s) = corr(rdm_roi1_HI(:),rdm_roi2_HI(:),'rows','pairwise','type','Spearman');
%             corrRoiRoi_LO(r1,r2,s) = corr(rdm_roi1_LO(:),rdm_roi2_LO(:),'rows','pairwise','type','Spearman');
            corrRoiRoi_HI(r1,r2,s) = rankCorr_Kendall_taua(rdm_roi1_HI(:),rdm_roi2_HI(:));
            corrRoiRoi_LO(r1,r2,s) = rankCorr_Kendall_taua(rdm_roi1_LO(:),rdm_roi2_LO(:));
        end
    end
end, clear r1 r2

figure('color',[1 1 1])
roiNamesTrue = {'Lingual','lHPC','rHPC','ACC','rIns','vmPFC','pOFC','aOFC'};

% plot val HI
corrRoiRoi_HI_mean = mean(corrRoiRoi_HI,3);
subplot(2,2,1),imagesc(corrRoiRoi_HI_mean,[0.07 0.18])
set(gca,'XTick',1:numel(roiNames),'fontsize',14,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',14,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45),title(['High ',whatScore])

% plot val LO
corrRoiRoi_LO_mean = mean(corrRoiRoi_LO,3);
subplot(2,2,2),imagesc(corrRoiRoi_LO_mean,[0.07 0.18])
set(gca,'XTick',1:numel(roiNames),'fontsize',14,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',14,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45),title(['Low ',whatScore])

% plot difference
subplot(2,2,3),imagesc(corrRoiRoi_HI_mean - corrRoiRoi_LO_mean)
set(gca,'XTick',1:numel(roiNames),'fontsize',14,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',14,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45),title('Difference')

% ttest
for i = 1:numel(roiNames)
    for j = 1:numel(roiNames)
        roiHI = corrRoiRoi_HI(i,j,:);
        roiLO = corrRoiRoi_LO(i,j,:);
        [~,p(i,j),ci,stats] = ttest(roiHI-roiLO,0,'tail','right');
        tfoo = stats.tstat;
        tvalues(i,j) = tfoo;
    end
end, clear roiLO roiHI tfoo stats corrRoiRoi_LO corrRoiRoi_HI

% plot difference
foo = p < 0.05;
pSignificant = tvalues.*foo;
subplot(2,2,4),imagesc(pSignificant,[-3 3])
set(gca,'XTick',1:numel(roiNames),'fontsize',14,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',14,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45),title('p < 0.05')
clear foo pSignificant roiNamesTrue


%% dimensionality
dimThreshold = 75;
for r = 1:numel(roiNames)
    for s = 1:numel(subs)
        explained = [];
        
        % subject name
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        dataMatrix = B_struct.(roiNames{r}).(subjName)';
        
        % reorder according to session
        %         newOrd = ordData(subs(s)).norm2sessions;
        %         dataMatrix = dataMatrix(newOrd,:);
        
        % take values in this session
        valuesSession = bData(subs(s)).imagination.val;
        valuesSession = valuesSession(:) + 0.00001*[1:60]';
        valuesSession_LO = valuesSession < median(valuesSession);
        valuesSession_HI = valuesSession > median(valuesSession);
        
        [~,~,~,~,explained(:,1),~] = pca(dataMatrix(1:60,:));
        [~,~,~,~,explained(:,2),~] = pca(dataMatrix(61:120,:));
        [~,~,~,~,explained(:,3),~] = pca(dataMatrix(121:180,:));
        [~,~,~,~,explained(:,4),~] = pca(dataMatrix(181:240,:));
        
        foo = cumsum(explained);
        
        for i = 1:4
            dimSession.(roiNames{r})(s,i) = find(foo(:,i) > dimThreshold,1);
        end
    end
    dimSession_mean.(roiNames{r}) = mean(dimSession.(roiNames{r}),1);
end, clear subjName i r s

%% Compare good and bad subjects
% subsBest = [15,2,23,16,1,25,29,18,21,11,12,22,4,6,19,24];
% subsWors = [30,28,9,14,20,17,13,27,28,11,10,3,8,5,7];
%
% corrRoiRoi_mean_best = mean(corrRoiRoi(:,:,subsBest),3);
% figure,imagesc(corrRoiRoi_mean_best,[0 0.18])
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
%
% corrRoiRoi_mean_wors = mean(corrRoiRoi(:,:,subsWors),3);
% figure,imagesc(corrRoiRoi_mean_wors,[0 0.18])
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
%
% corrRoiRoi_mean_diff = corrRoiRoi_mean_best - corrRoiRoi_mean_wors;
% figure,imagesc(corrRoiRoi_mean_diff)
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
%
% for i = 1:numel(roiNames)
%     for j = 1:numel(roiNames)
%         bbb = corrRoiRoi(i,j,subsBest);
%         www = corrRoiRoi(i,j,subsWors);
%         [~,p(i,j),ci,stats] = ttest(bbb-www);
%         ttt = stats.tstat;
%         aaa(i,j) = ttt;
%     end
% end
% figure,imagesc(aaa)
% set(gca,'XTick',1:numel(roiNames),'fontsize',11,'XtickLabel',roiNamesTrue,...
%     'YTick',1:numel(roiNames),'fontsize',11,'YtickLabel',roiNamesTrue)
% xtickangle(45),ytickangle(45)
