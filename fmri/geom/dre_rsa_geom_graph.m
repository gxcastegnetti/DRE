%% dre_rsa_geom_graph
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_pulse_FvB';

%% directories
fs         = filesep;
dir.geoCod = pwd;
idcs       = strfind(dir.geoCod,'/');
dir.dre    = dir.geoCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsaCod = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];

% paths
addpath([pwd,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 50];

%% properties
propNames = {'distAvg','distVar','distSke','distKur','distVarCoeff'};
propNames = {'corrFB'};

%% loop over subjects
if true
    for s = 1:length(subs)
        
        % load metadata from SPM
        dirFun = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
        d = spm_select('List', dirFun, '^uaf.*\.nii$');
        d = d(end-1,:);
        epi_file = {[dirFun fs d]};
        subjectMetadataStruct = spm_vol(epi_file);
        
        % load geometry maps
        load([dir.out,fs,analysisName,fs,'sl_SF',num2str(subs(s),'%03d')],'geomDiff');
        
        %% compute properties
        %         distProp{1} = mean(geomDiff,4);
        %         distProp{2} = var(geomDiff,0,4);
        %         distProp{3} = skewness(geomDiff,1,4);
        %         distProp{4} = kurtosis(geomDiff,1,4);
        %         distProp{5} = distProp{2}./distProp{1};
        distProp{1} = geomDiff;
        
        %% loop over properties
        for p = 1:length(propNames)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Write the native-space r-map to a file %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            gMapMetadataStruct_nS = subjectMetadataStruct{1};
            gMapMetadataStruct_nS.fname = [dir.out,fs,analysisName,fs,propNames{p},fs,'Map_',propNames{p},'_SF',num2str(subs(s),'%03d'),'.nii'];
            gMapMetadataStruct_nS.descrip =  propNames{p};
            gMapMetadataStruct_nS.dim = size(distProp{p});
            gMapMetadataStruct_nS.dt = [16 0];
            gotoDir([dir.out,fs,analysisName,fs,propNames{p}]);
            spm_write_vol(gMapMetadataStruct_nS, distProp{p});
            
            %%%%%%%%%%%%%%%%%%%%
            % normalise to MNI %
            %%%%%%%%%%%%%%%%%%%%
            
            % select forward deformation images from T1 segmentation step
            dirStruct = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
            d = spm_select('List', dirStruct, '^y_.*\.nii$');
            y_file = {[dirStruct fs d]};
            
            % ------- g-map -------
            
            % select correlation map
            gMapFile = {gMapMetadataStruct_nS.fname};
            
            % prepare job for correlation map
            job_rs{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            job_rs{1}.spatial{1}.normalise{1}.write.subj.resample = gMapFile;
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job_rs{1}.spatial{1}.normalise{1}.write.woptions.interp = 1; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d'),' - g-map ',propNames{p}])
            d = spm_jobman('run',job_rs);
            clear job
            
            % read them back in
            wMaskFile = [dir.rsaOut,fs,'sl',fs,'_nS_masks_gm',fs,'wnS_gm_SF',num2str(subs(s),'%03d'),'.nii'];
            mask_sS = spm_read_vols(spm_vol(wMaskFile));
            
            %%%%%%%%%%
            % smooth %
            %%%%%%%%%%
            
            % unsmoothed file
            wrMapFile = [dir.out,fs,analysisName,fs,propNames{p},fs,'wMap_',propNames{p},'_SF',num2str(subs(s),'%03d'),'.nii'];
            
            disp(['Smoothing sub#', num2str(subs(s),'%03d'),' - ',propNames{p}])
            swMapFile = [dir.out,fs,analysisName,fs,propNames{p},fs,'swMap_',propNames{p},'_SF',num2str(subs(s),'%03d'),'.nii'];
            spm_smooth(wrMapFile,swMapFile,[3 3 3]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % re-apply mask after smoothing %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % take smoothed-norm and norm data
            swrMetadataStruct_sS = spm_vol(swMapFile);
            wrMetadataStruct_sS = spm_vol(wrMapFile);
            swrData = spm_read_vols(swrMetadataStruct_sS);
            wrData = spm_read_vols(wrMetadataStruct_sS);
            
            % apply mask
            swrData(mask_sS == 0) = nan;
            wrData(mask_sS == 0) = nan;
            
            % write again
            spm_write_vol(swrMetadataStruct_sS, swrData);
            spm_write_vol(wrMetadataStruct_sS, wrData);
            
        end
    end
end

%% load and concatenate files
for s = 1:length(subs)
    
    % loop over models
    for p = 1:length(propNames)
        
        % select file
        swMapFile = [dir.out,fs,analysisName,fs,propNames{p},fs,'wMap_',propNames{p},'_SF',num2str(subs(s),'%03d'),'.nii'];
        
        % mask smoothed r-maps with MNI mask
        swMap = spm_read_vols(spm_vol(swMapFile));
        
        % concatenate across subjects
        rMaps_all.(propNames{p})(:,:,:,s) = swMap;
        
        %         for i = 1:79
        %             figure
        %             subplot(1,2,1),imagesc(swMap(:,:,i))
        %         end
        
    end
end

%% statistics
for p = 1:length(propNames)
    
    %%%%%%%%%%%%%%%%%
    % compute g-map %
    %%%%%%%%%%%%%%%%%
    
    % take r-map for current model
    gMaps = rMaps_all.(propNames{p});
    
    impMask = ~isnan(sum(gMaps,4));
    
    % obtain a pMaps from applying a 1-sided signrank test and also t-test to
    % the model similarities:
    p1 = NaN(79,95,79);
    p2 = p1;
    t1 = p1;
    for x = 1:size(gMaps,1)
        for y = 1:size(gMaps,2)
            for z = 1:size(gMaps,3)
                foo = squeeze(gMaps(x,y,z,:));
                foo(isnan(foo)) = [];
                if length(foo) < 15, continue, end
                [~, p1(x,y,z), ~, stats] = ttest(squeeze(gMaps(x,y,z,:)));
                t1(x,y,z) = stats.tstat;
%                 [p2(x,y,z)] = signrank_onesided(squeeze(gMaps(x,y,z,:)));
                
            end
        end
        disp(x);
    end
    
    % update user
    disp(['property: ',propNames{p}])
    
    %     % apply FDR correction
    %     pThrsh_t  = FDRthreshold(p1,0.05,impMask);
    %     pThrsh_sr = FDRthreshold(p2,0.05,impMask);
    
    % mark the suprathreshold voxels in yellow
    %     supraThreshMarked_t = zeros(size(p1));
    %     supraThreshMarked_t(p1 <= pThrsh_t) = 1;
    %     supraThreshMarked_sr = zeros(size(p2));
    %     supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;
    
    % write g-map
    swMapFile = [dir.out,fs,analysisName,fs,propNames{p},fs,'swMap_',propNames{p},'_SF',num2str(subs(s),'%03d'),'.nii'];
    tMapMetadataStruct_sS = spm_vol(swMapFile);
    tMapMetadataStruct_sS.fname = [dir.out,fs,analysisName,fs,propNames{p},fs,'neg_meanMap_',propNames{p},'.nii'];
    tMapMetadataStruct_sS.descrip = 'g-map';
    tMapMetadataStruct_sS.dim = size(t1);
    spm_write_vol(tMapMetadataStruct_sS, mean(gMaps,4));
    
    % write t-map
    tMapMetadataStruct_sS = spm_vol(swMapFile);
    tMapMetadataStruct_sS.fname = [dir.out,fs,analysisName,fs,propNames{p},fs,'neg_tMap_',propNames{p},'.nii'];
    tMapMetadataStruct_sS.descrip = 't-map';
    tMapMetadataStruct_sS.dim = size(t1);
    spm_write_vol(tMapMetadataStruct_sS, t1);
    
    % write p-map
    %     pMapMetadataStruct_sS = spm_vol(swrMapFile);
    %     pMapMetadataStruct_sS.fname = [dirSl,fs,modelName,fs,'pMap_',modelName,'.nii'];
    %     pMapMetadataStruct_sS.descrip = 'p-map';
    %     pMapMetadataStruct_sS.dim = size(supraThreshMarked_t);
    %     supraThreshMarked_t = p1 < 0.005;
    %     spm_write_vol(pMapMetadataStruct_sS, supraThreshMarked_t);
    
end
