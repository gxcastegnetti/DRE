%% dre_rsa_sl_test
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
% analysisName = 'rsa_sl_pulse_ons0';
analysisName = 'rsa_sl_pulse_ons0';
% analysisName = 'dim_sl_ons0';
% betaid       = 'rsa_pulse_ons0';
betaid       = 'rsa_pulse_ons0';
thisIsDim    = false;

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.dimOut = [dir.dre,fs,'out',fs,'fmri',fs,'dim'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];

% paths
dir.spm  = '/Users/gcastegnetti/Desktop/tools/matlab/spm12';
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath(dir.spm))

%% Subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
% subsBest = sort([23 18 5 3 21 11 10 20 17 28 24  1 15 22]);
% subsWors = sort([2  14 6 4  7  9 27 26 12 16 19 13  8 25]);
% subs = subs(subsBest);

%% Set options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% mask <- should be 79x95x79 because we put correlation map in MNI

% reslice mask to the size of the correlation maps, if needed
if ~exist([dir.mskOut,fs,'rgm.nii'],'file')
    
    % load sample EPI from one subject for coregistration
    dirFun = [dir.datScn,fs,'SF039',fs,'fun',fs,'S4'];
    d = spm_select('List', dirFun, '^wuaf.*\.nii$');
    d = d(end-1,:);
    epi_file = [dirFun fs d];
    job_rs{1}.spm.spatial.coreg.write.ref = {epi_file};
    
    % select mask to coregister
    job_rs{1}.spm.spatial.coreg.write.source = {[dir.mskOut,fs,'gm.nii']};
    
    % defaults
    job_rs{1}.spm.spatial.coreg.write.roptions.interp = 4;
    job_rs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    job_rs{1}.spm.spatial.coreg.write.roptions.mask = 0;
    job_rs{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    % run job
    disp('Coregistering grey matter mask')
    spm_jobman('run',job_rs);
    clear job d epi_file
end

mask = niftiread([dir.mskOut,fs,'rgm.nii']);
mask = logical(mask);

%% model names
modelNames = {'val','fam','oid','cxt'};
% modelNames = {'dval','vCho','vUnc','chMu','oids'};
% modelNames = {'valMed','famMed'};

if thisIsDim
    modelNames = {'dim'};
end
%% soecify some directories

% directory with searchlight correlation maps
dirSl = [userOptions.rootPath,filesep,'sl',fs,analysisName];

if thisIsDim
    % directory with searchlight correlation maps
    dirSl = [dir.dimOut,fs,'sl',fs,analysisName];
end

% directory with betas
dirBeta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];

%% loop over subjects
if true
    for s = 1:length(subs)
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % make r-maps SPM-like %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % load metadata from SPM
        dirFun = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
        d = spm_select('List', dirFun, '^uaf.*\.nii$');
        d = d(end-1,:);
        epi_file = {[dirFun fs d]};
        betaFile = [dirBeta,fs,'SF',num2str(subs(s),'%03d'),fs,'beta_0001.nii'];
        subjectMetadataStruct = spm_vol(epi_file);
        
        % load correlation maps
        %     load([dirSl,fs,'sl_context_SF',num2str(subs(s),'%03d'),'.mat']);
        if ~thisIsDim
            load([dirSl,fs,'sl_SF',num2str(subs(s),'%03d'),'.mat']);
        else
            load([dirSl,fs,'sl_dim_SF',num2str(subs(s),'%03d'),'.mat']);
        end
        
        %% loop over models
        for m = 1:length(modelNames)
            
            % current model name
            modelName = modelNames{m};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Write the native-space r-map to a file %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            rMapMetadataStruct_nS = subjectMetadataStruct{1};
            rMapMetadataStruct_nS.fname = [dirSl,fs,modelName,fs,'rMap_',modelName,'_SF',num2str(subs(s),'%03d'),'.nii'];
            rMapMetadataStruct_nS.descrip =  'R-map';
            rMapMetadataStruct_nS.dim = size(rs(:,:,:,m));
            rMapMetadataStruct_nS.dt = [16 0];
            gotoDir([dirSl,fs,modelName]);
            spm_write_vol(rMapMetadataStruct_nS, rs(:,:,:,m));
            
            %%%%%%%%%%%%%%%%%%%%
            % normalise to MNI %
            %%%%%%%%%%%%%%%%%%%%
            
            % select forward deformation images from T1 segmentation step
            dirStruct = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
            d = spm_select('List', dirStruct, '^y_.*\.nii$');
            y_file = {[dirStruct fs d]};
            
            % ------- r-map -------
            
            % select correlation map
            rs_file = {rMapMetadataStruct_nS.fname};
            
            % prepare job for correlation map
            job_rs{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            job_rs{1}.spatial{1}.normalise{1}.write.subj.resample = rs_file;
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.interp = 1; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d'),' - R-map ',modelName])
            d = spm_jobman('run',job_rs);
            clear job
            
            % ------- mask -------
            
            % write the native-space mask to a file
            maskMetadataStruct_nS = subjectMetadataStruct{1};
            maskMetadataStruct_nS.fname = [dir.out,fs,'sl',fs,'_nS_masks_gm',fs,'nS_gm_SF',num2str(subs(s),'%03d'),'.nii'];
            maskMetadataStruct_nS.descrip =  'Native space mask';
            maskMetadataStruct_nS.dim = size(mask);
            spm_write_vol(maskMetadataStruct_nS, mask);
            
            % select mask file we just saved
            mask_file = {maskMetadataStruct_nS.fname};
            
            % prepare job for mask
            job_mask{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            job_mask{1}.spatial{1}.normalise{1}.write.subj.resample = mask_file;
            job_mask{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job_mask{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job_mask{1}.spatial{1}.normalise{1}.write.woptions.interp = 1; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d'),' - MASK ',modelName])
            d = spm_jobman('run',job_mask);
            clear job
            
            % read them back in
            wMaskFile = [dir.out,fs,'sl',fs,'_nS_masks_gm',fs,'wnS_gm_SF',num2str(subs(s),'%03d'),'.nii'];
            wrMapFile = [dirSl,fs,modelName,fs,'wrMap_',modelName,'_SF',num2str(subs(s),'%03d'),'.nii'];
            mask_sS = spm_read_vols(spm_vol(wMaskFile));
            
            % Fix the normalisation of the mask
            maskMetadataStruct_sS = spm_vol(wMaskFile);
            maskMetadataStruct_sS.fname = wMaskFile;
            maskMetadataStruct_sS.descrip =  'Common space mask';
            maskThreshold = 0.01;
            mask_sS(mask_sS < maskThreshold) = 0;
            mask_sS(isnan(mask_sS)) = 0;
            maskMetadataStruct_sS.dim = size(mask_sS);
            spm_write_vol(maskMetadataStruct_sS, mask_sS);
            
            %%%%%%%%%%
            % smooth %
            %%%%%%%%%%
            
            if ~thisIsDim
                disp(['Smoothing sub#', num2str(subs(s),'%03d'),' - ',modelName])
                swrMapFile = [dirSl,fs,modelName,fs,'swrMap_',modelName,'_SF',num2str(subs(s),'%03d'),'.nii'];
                
                spm_smooth(wrMapFile,swrMapFile,[9 9 9]);
            end
            
        end
    end
end

%% load and concatenate files
for s = 1:length(subs)
    
    % loop over models
    for m = 1:length(modelNames)
        
        % current model name
        modelName = modelNames{m};
        
        % select file
        swrMapFile = [dirSl,fs,modelName,fs,'swrMap_',modelName,'_SF',num2str(subs(s),'%03d'),'.nii'];
        
        % mask smoothed r-maps with MNI mask
        swrMap = spm_read_vols(spm_vol(swrMapFile));
        swrMap(mask == 0) = NaN;
        
        % concatenate across subjects
        rMaps_all.(modelName)(:,:,:,s) = swrMap;
        % figure,imagesc(swrMap(:,:,40)),colorbar
        
    end
end

%% if it's a dimensionality searchlight and plot
if thisIsDim
    meanDim = mean(rMaps_all.dim,4);
    meanDim(meanDim <= 1) = nan;
    
    %     for i = 1:79
    %         figure,imagesc(meanDim(:,:,i)),colorbar
    %     end
    
    % write d(imensionality)-map
    dMapMetadataStruct_sS = spm_vol(swrMapFile);
    dMapMetadataStruct_sS.fname = [dirSl,fs,'dMap.nii'];
    dMapMetadataStruct_sS.descrip =  'd-map';
    dMapMetadataStruct_sS.dim = size(meanDim);
    spm_write_vol(dMapMetadataStruct_sS, meanDim);
    
    % end code here
    return
end

%% statistics
for m = 1:length(modelNames)
    
    modelName = modelNames{m};
    
    %%%%%%%%%%%%%%%%%
    % compute p-map %
    %%%%%%%%%%%%%%%%%
    
    % take r-map for current model
    rMaps = rMaps_all.(modelName);
    
    % obtain a pMaps from applying a 1-sided signrank test and also t-test to
    % the model similarities:
    p1 = NaN(79,95,79);
    p2 = p1;
    t1 = p1;
    for x = 1:size(rMaps,1)
        for y = 1:size(rMaps,2)
            for z = 1:size(rMaps,3)
                if mask(x,y,z) == 1
                    [~, p1(x,y,z), ~, stats] = ttest(squeeze(rMaps(x,y,z,:)),0,0.05,'right');
                    if p1(x,y,z) < 0.005
                        t1(x,y,z) = stats.tstat;
                    end
                    [p2(x,y,z)] = signrank_onesided(squeeze(rMaps(x,y,z,:)));
                end
            end
        end
        disp(x);
    end
    
    % plot p-values
    if true && m == 25
        for i = 1:79
            figure,
            subplot(1,2,1),imagesc(p1(:,:,i))
            subplot(1,2,2),imagesc(p1(:,:,i)<0.05)
        end
    end
    
    % apply FDR correction
    pThrsh_t  = FDRthreshold(p1,0.05,mask);
    pThrsh_sr = FDRthreshold(p2,0.05,mask);
    
    % mark the suprathreshold voxels in yellow
    supraThreshMarked_t = zeros(size(p1));
    supraThreshMarked_t(p1 <= pThrsh_t) = 1;
    supraThreshMarked_sr = zeros(size(p2));
    supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;
    
    % write t-map
    swrMapFile = [dirSl,fs,modelName,fs,'swrMap_',modelName,'_SF',num2str(subs(s),'%03d'),'.nii'];
    tMapMetadataStruct_sS = spm_vol(swrMapFile);
    tMapMetadataStruct_sS.fname = [dirSl,fs,modelName,fs,'tMap_',modelName,'.nii'];
    tMapMetadataStruct_sS.descrip = 't-map';
    tMapMetadataStruct_sS.dim = size(supraThreshMarked_t);
    spm_write_vol(tMapMetadataStruct_sS, t1);
    
    % write p-map
    pMapMetadataStruct_sS = spm_vol(swrMapFile);
    pMapMetadataStruct_sS.fname = [dirSl,fs,modelName,fs,'pMap_',modelName,'.nii'];
    pMapMetadataStruct_sS.descrip = 'p-map';
    pMapMetadataStruct_sS.dim = size(supraThreshMarked_t);
    supraThreshMarked_t = p1 < 0.005;
    spm_write_vol(pMapMetadataStruct_sS, supraThreshMarked_t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create cluster-specific masks %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sum(supraThreshMarked_sr(:)) > 0 && pThrsh_sr < Inf && false % <------- false
        
        % find connected clusters
        k = 1;
        cluster = {};
        L = bwlabeln(supraThreshMarked_sr);
        
        % loop over such clusters
        for i = 1:max(L(:))
            cluster_foo{i} = L == i; %#ok<*SAGROW>
            clusterSize(i) = sum(cluster_foo{i}(:));
            
            % take only clusters with size > number of conditions per session (60)
            if clusterSize(i) >= 60
                cluster{k} = cluster_foo{i};
                
                % create and save mask
                maskMetadataStruct_sS = pMapMetadataStruct_sS;
                maskMetadataStruct_sS.fname = [dirSl,fs,'mask_sl',modelName,'_',num2str(k),'.nii'];
                spm_write_vol(maskMetadataStruct_sS, cluster{k});
                
                % update index
                k = k+1;
            end
        end, clear cluster_foo i k L
    end
    
end


