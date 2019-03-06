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
roiNames = {'rsaVal_LG_10mm','rsaVal_ITG','rsaVal_PCC_10mm','l_hpc','r_hpc','rsaVal_ACC_10mm','rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm'};
roiNames = {'rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm','rsaCon_vmPFC_10mm'};
roiNames = {'l_hpc','r_hpc','rsaVal_vmPFC_10mm'};
% roiNamesTrue = {'LG','ITG','PCC','l HPC','r HPC','ACC','vmPFC','OFC','dlPFC'};
roiNamesTrue = {'lHPC','rHPC','vmPFC'};

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
        foo = foo + RDM_brain(r,s).RDM/numel(subs);
    end
    
    % fill struct for the RSA toolbox
    RDM_MDS{r}.color = [0 0 1];
    RDM_MDS{r}.name = roiNames{r};
    RDM_MDS{r}.RDM = foo;
    
    % plot
    %     subplot(3,5,r)
    %     RDM_MDS(r).RDM(RDM_MDS(r).RDM==0)=nan;
    %     foo = nanmean(RDM_MDS(r).RDM(:));
    %     imagesc(RDM_MDS(r).RDM),title(roiNames{r}),caxis([foo-1 foo+1])
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
fooVal = zeros(size(B,1),size(B,1));
fooCon = zeros(size(B,1),size(B,1));
if size(RDM_brain(1,1).RDM,1) == 240
    for s = 1:length(subs)
        RDM_model(1,s).name = 'val';
        RDM_model(1,s).RDM = RDMs_models{s}.val; %#ok<*SAGROW>
%         RDM_model(1,s).color = [0 1 0];
%         RDM_model(2,s).name = 'con';
%         RDM_model(2,s).RDM = RDMs_models{s}.con;
        RDM_model(2,s).color = [0 1 0];
        RDM_model(2,s).name = 'oid';
        RDM_model(2,s).RDM = 1-mat_ID;
%         RDM_model(3,s).color = [0 1 0];
%         RDM_model(3,s).name = 'oid';
%         RDM_model(3,s).RDM = 50*RDMs_models{s}.con+RDMs_models{s}.val;
%         RDM_model(3,s).color = [0 1 0];
        
        % averages
        fooVal = fooVal + RDM_model(1,s).RDM/numel(subs);
        fooCon = fooCon + RDM_model(2,s).RDM/numel(subs);
    end
    
    % for every region and sub, correlate RDM and model
    scoreNames = {'Value','Confidence','Familiarity','Object ID','Context'};
    for r = 1:length(roiNames)
        for s = 1:size(RDM_brain,2)
            for m = 1:size(RDM_model,1)
                a = vectorizeRDM(RDM_brain(r,s).RDM);
                b = vectorizeRDM(RDM_model(m,s).RDM);
                corrRoiModel(r,s,m) = corr(a(:),b(:),'rows','pairwise','type','spearman');
                %                 rL2(m,s) = fisherTransform(rL2(m,s));
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

%%%%%%%%%%%%
% MDS plot %
%%%%%%%%%%%%

% RDM_MDS{numel(roiNames)+1}.color = [0 0 1];
% RDM_MDS{numel(roiNames)+1}.name = 'Value';
% RDM_MDS{numel(roiNames)+1}.RDM = fooVal;
% RDM_MDS{numel(roiNames)+2}.color = [0 0 1];
% RDM_MDS{numel(roiNames)+2}.name = 'Confidence';
% RDM_MDS{numel(roiNames)+2}.RDM = fooCon;
% userOptions.rubberbands = true;
% MDSRDMs(RDM_MDS,userOptions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ttest of correlations for each ROI and model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(roiNames)
    for m = 1:size(RDM_model,1)
        scores = corrRoiModel(r,:,m);
        %         [h,pCorr(r,m),~,~] = ttest(scores,0,'Tail','right');
        [pCorr(r,m),h,~] = signrank(scores,0,'Tail','right');
        [rPerf.r(r,m),rPerf.p(r,m)] = corr(scores',slopes','type','spearman');
    end
end, clear r m mat_ID

%%%%%%%%%%%%%%%%%%%%%
% plot mean and sem %
%%%%%%%%%%%%%%%%%%%%%
means = squeeze(mean(corrRoiModel(:,:,1:2),2));
sems  = squeeze(std(corrRoiModel(:,:,1:2),0,2)/sqrt(numel(subs)));

% means = means(end-1:end,1:2);
% sems = sems(end-1:end,1:2);

roiNamesTrue = {'vmPFC','lOFC','dlPFC','mOFC'};
myColors = [0,128,255; 255,51,153]/255;
myColors_ss = [170,200,255; 255,204,229]/255;
figure('color',[1 1 1])

hb = bar(means,0.6); hold on
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
        
    for jb = 1:numel(roiNames)
        scatter(xData(jb)+0.03*randn(numel(subs),1),corrRoiModel(jb,:,ib),25,'MarkerEdgeColor',myColors_ss(ib,:),...
            'MarkerFaceColor',myColors_ss(ib,:))
    end
    errorbar(xData',means(:,ib),sems(:,ib),'linestyle','none','color','k','linewidth',2.5,'capsize',0)
    hb(ib).FaceColor = myColors(ib,:);
end, clear ib

legend({'Value','Confid.'},'location','southeast'), legend boxoff
set(gca,'fontsize',18,'xtick',1:numel(roiNamesTrue),...
    'xticklabels',roiNamesTrue), ylim([-0.025 0.0301]),xtickangle(45)
ylabel('Correlation(ROI, model)')


%% ROI-ROI comparison in terms of correlation with model
for r1 = 1:length(roiNames)
    for r2 = 1:length(roiNames)
        corrRoiVal_r1 = squeeze(corrRoiModel(r1,:,1));
        corrRoiVal_r2 = squeeze(corrRoiModel(r2,:,1));
        %         [~,p(r1,r2),~,~] = ttest(corrRoiVal_r2 - corrRoiVal_r1,0,'tail','right');
        [pRoiRoi(r1,r2),~,~] = signrank(corrRoiVal_r2 - corrRoiVal_r1,0,'tail','right');
    end
end


%% ROI-ROI correlations

whatScore = 'confidence';

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
    pctile_con_33 = prctile(Y_con,25);
    pctile_con_66 = prctile(Y_con,75);
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
% roiNamesTrue = {'lHPC','rHPC','vmPFC','OFC','dlPFC'};

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
foo = p < 0.01;
pSignificant = tvalues.*foo;
subplot(2,2,4),imagesc(pSignificant,[-3 3])
set(gca,'XTick',1:numel(roiNames),'fontsize',14,'XtickLabel',roiNamesTrue,...
    'YTick',1:numel(roiNames),'fontsize',14,'YtickLabel',roiNamesTrue)
xtickangle(45),ytickangle(45),title('p < 0.01')
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
