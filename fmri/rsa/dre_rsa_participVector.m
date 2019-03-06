%% dre_rsa_participVector
% ~~~
% GX Castegnetti --- 2019

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
subs = [4 5 7:9 13:17 19 21 23 25:26 29:32 34 35 37 39 40:43 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,3)];

%% extract behavioural data and rearrange for visualisation
bData = dre_extractData(dir,subs,taskOrd,0);
ordData = dre_rearrange_3L(dir,subs,taskOrd,bData);

%% which mask?

roiNames = {'rsaVal_LG_10mm','rsaVal_ACC_10mm','rsaVal_vmPFC_10mm','rsaVal_OFC_10mm','rsaVal_dlPFC_10mm'};
roiNames = {'rsaVal_dlPFC_10mm'};

%% prewhiten activity in the mask
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        % subject name
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        % SPM file from 1st level analysis
        subjSPMFile = [dir.beta,fs,'SF',num2str(subs(s),'%03d'),fs,'SPM.mat'];
        load(subjSPMFile)
        
        % subjective mask
        subjMaskFile.fname = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        
        % load indices for normal ordering
        fooDir = [dir.dre,filesep,'out',filesep,'fmri',filesep,'rsa'];
        load([fooDir,filesep,'toNormalOrder',filesep,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
        
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
            B = B(toNormalOrder,:);
        end
        
        % construct RDM
        rdm = squareform(pdist(B,'correlation'));
        RDM_brain(r,s).RDM = rdm;
    end
end, clear rdm fooDir subjSPMFile s SPM betaid

%% FIRE vs BOAT
for s = 1:length(subs)
    
    disp(['sub#',num2str(subs(s))])
    
    % load indices for normal ordering
    fooDir = [dir.dre,filesep,'out',filesep,'fmri',filesep,'rsa'];
    load([fooDir,filesep,'toNormalOrder',filesep,'SF',num2str(subs(s),'%03d')],'toNormalOrder')
    
    % extract presentation indices
    objVal = [bData(subs(s)).imagination(1).val; bData(subs(s)).imagination(2).val;...
        bData(subs(s)).imagination(3).val; bData(subs(s)).imagination(4).val];
    objCon = [bData(subs(s)).imagination(1).con; bData(subs(s)).imagination(2).con;...
        bData(subs(s)).imagination(3).con; bData(subs(s)).imagination(4).con];
    
    % sort indices
    objVal = objVal(toNormalOrder);
    objCon = objCon(toNormalOrder);
    %     objVal_F = objVal(1:120);
    %     objVal_B = objVal(121:end);
    
    % brain activities during fire and boat trials
    %     X_F = B(:,1:120)';
    %     X_B = B(:,121:end)';
    %
    %     % remove nans
    %     objVal_F_real = objVal_F(~isnan(objVal_F));
    %     objVal_B_real = objVal_B(~isnan(objVal_B));
    %
    %     % loop over voxels
    %     numVoxels = size(B,1);
    %     for vx = 1:numVoxels
    %
    %         % fire
    %         y_F = objVal_F_real;
    %         x_F = X_F(vx,~isnan(objVal_F));
    %         XX_F = [ones(numel(x_F),1) x_F(:)];
    %         betas_F = XX_F\y_F;
    %         slopes_F(vx) = abs(betas_F(2));
    %
    %         % boat
    %         y_B = objVal_B_real;
    %         x_B = X_B(vx,~isnan(objVal_B));
    %         XX_B = [ones(numel(x_B),1) x_B(:)];
    %         betas_B = XX_B\y_B;
    %         slopes_B(vx) = abs(betas_B(2));
    %     end
    %     r_betas(s) = corr(slopes_F(:),slopes_B(:),'type','spearman','rows','complete');
    %
    %     % permutation test
    %     nPerm = 100;
    %     for i = 1:nPerm
    %
    %         % permute values
    %         fooRand_val_F = randperm(numel(objVal_F));
    %         fooRand_val_B = randperm(numel(objVal_B));
    %         objVal_F_perm = objVal(fooRand_val_F);
    %         objVal_B_perm = objVal(fooRand_val_B);
    %
    %         % remove nans
    %         objVal_F_perm = ~isnan(objVal_F_perm);
    %         objVal_B_perm = ~isnan(objVal_B_perm);
    %
    %         for vx = 1:numVoxels
    %
    %             % fire
    %             y_F_perm = objVal_F_perm;
    %             x_F = X_F(vx,:);
    %             XX_F = [ones(numel(x_F),1) x_F(:)];
    %             betas_F_perm = XX_F\y_F_perm;
    %             slopes_F_perm(vx) = abs(betas_F_perm(2));
    %
    %             % boat
    %             y_B_perm = objVal_B_perm;
    %             x_B = X_B(vx,:);
    %             XX_B = [ones(numel(x_B),1) x_B(:)];
    %             betas_B_perm = XX_B\y_B_perm;
    %             slopes_B_perm(vx) = abs(betas_B_perm(2));
    %
    %         end
    %
    %         r_betas_perm(s,i) = corr(slopes_F_perm(:),slopes_B_perm(:),'type','spearman','rows','complete');
    %         [pVal_perm(i),~,~] = signrank(r_betas_perm(:,i));
    %     end
    
    %% value/confidence
    numVoxels = size(B,2);
    for vx = 1:numVoxels
        
        % voxel activity
        y = B(:,vx);
        
        % value
        X_val = [ones(length(objVal),1) objVal];
        betas_val = X_val(~isnan(objVal),:)\y(~isnan(objVal));
        slopes_val(vx) = abs(betas_val(2));
        
        % confidence
        X_con = [ones(length(objCon),1) objCon];
        betas_con = X_con(~isnan(objCon),:)\y(~isnan(objCon));
        slopes_con(vx) = abs(betas_con(2));
        
    end
    r_betas(s) = corr(slopes_val',slopes_con','type','spearman');
    
    % permutation test
    nPerm = 1000;
    for i = 1:nPerm
        fooRand_val = randperm(numel(objVal));
        fooRand_con = randperm(numel(objCon));
        objVal_perm = objVal(fooRand_val);
        objCon_perm = objCon(fooRand_con);
        
        for vx = 1:numVoxels
            
            % voxel activity
            y = B(:,vx);
            
            % value
            X_val = [ones(length(objVal_perm),1) objVal_perm];
            betas_val = X_val(~isnan(objVal_perm),:)\y(~isnan(objVal_perm));
            slopes_val_perm(vx) = abs(betas_val(2));
            
            % confidence
            X_con = [ones(length(objCon_perm),1) objCon_perm];
            betas_con = X_con(~isnan(objCon_perm),:)\y(~isnan(objCon_perm));
            slopes_con_perm(vx) = abs(betas_con(2));
            
        end
        
        r_betas_perm(s,i) = corr(slopes_val_perm',slopes_con_perm','type','spearman');
        %         [pVal_perm(i),~,~] = signrank(r_betas_perm(:,i));
    end, clear fooRand_val fooRand_con objVal_perm objCon_perm
    
end

% stats with real indices
[pVal_real,~,stats] = signrank(r_betas,0);
zVal_real = stats.zval;

for i = 1:nPerm
    [~,~,stats] = signrank(r_betas_perm(:,i),0);
    zVal_perm(i) = stats.zval;
end
zVal_larger = sum((zVal_perm) >= (zVal_real));
pVal = zVal_larger/nPerm;

%% plot
allSubs_lOFC = [0.162369597615499,-0.0782753014496681,0.188680395610351,0.0902519983742040,...
    -0.0242446822923723,0.0821230185611706,0.00502641918439236,0.109761549925484,-0.0194350359029942,...
    0.118242785530416,0.138111367023439,0.123323397913562,0.0119292778756266,-0.121955019645035,...
    0.0795014225714673,0.163798943232624,0.0159395745833898,0.0391681343991329,0.108068012464436,...
    -0.0229237230727544,0.0649708711556700,-0.00453868039561035,0.0699227746917762,0.0630266901503861,...
    0.0528248204850291,0.105107708982523,0.208061238314591,-0.0812288307817369,-0.0927652079664002,0.0282414307004471];

allSubs_mOFC = [-0.00126030951719025,-0.0208599759058475,0.0631266796404411,0.125549068668335,...
    0.0869057547956631,-0.0241682883884719,-0.0493096098600686,0.0953850430914651,-0.0412102678157724,...
    0.0613937540543045,0.0808729496802891,0.0285701047168937,0.0997034565841905,0.179696042998795,...
    -0.00979519970345659,-0.0732925586136595,-0.0682142526179223,0.0132610508757298,0.216791770920211,...
    -0.124205356315448,-0.0385043091465110,0.104550088036327,-0.117375590770086,-0.00807154109906404,...
    -0.0382355666759337,0.0292558613659531,0.219701603187842,-0.0473357427485868,-0.163237883421370,0.0258826800111204];

means(1) = mean(allSubs_lOFC);
means(2) = mean(allSubs_mOFC);

sems(1) = std(allSubs_lOFC)/sqrt(30);
sems(2) = std(allSubs_mOFC)/sqrt(30);

roiNamesTrue = {'lOFC','mOFC'};
myColors = [95,95,95]/255;
myColors_ss = [170,170,170]/255;
figure('color',[1 1 1])

hb = bar(means,0.4); hold on
hb.FaceColor = myColors;

scatter(1+0.03*randn(numel(subs),1),allSubs_lOFC,25,'MarkerEdgeColor',myColors_ss,...
    'MarkerFaceColor',myColors_ss)
scatter(2+0.03*randn(numel(subs),1),allSubs_mOFC,25,'MarkerEdgeColor',myColors_ss,...
    'MarkerFaceColor',myColors_ss)

errorbar(1,means(1),sems(1),'linestyle','none','color','k','linewidth',2.5,'capsize',0)
errorbar(2,means(2),sems(2),'linestyle','none','color','k','linewidth',2.5,'capsize',0)

set(gca,'fontsize',18,'xtick',1:2,'xticklabels',{'lOFC','mOFC'}),...
    ylim([-0.25 0.301]),xtickangle(45)
ylabel('Correlation(ROI, model)')

