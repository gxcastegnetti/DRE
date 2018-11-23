%% dre_rsa_roi_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_sphere';
% analysisName = 'rsa_roi_pulse_cxt';
dirBeta = 'rsa_pulse_ons0';

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox',fs,'rsa']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% load masks
% roiNames = {'sphere_9-0_-61_25','sphere_9-0_-43_35','sphere_9-0_-25_39','sphere_9-0_-7_40',...
%     'sphere_9-0_11_37','sphere_9-0_28_29','sphere_9-0_41_20','sphere_9-0_47_11'};

roiNames = {'box_w-16_16_16-0_-60_26','box_w-16_16_16-0_-44_36','box_w-16_16_16-0_-28_40','box_w-16_16_16-0_-12_42',...
    'box_w-16_16_16-0_4_42','box_w-16_16_16-0_20_36','box_w-16_16_16-0_36_23','box_w-16_16_16-0_46_7'};

roiNames = {'hpc_lr'};

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];

%% reverse normalise mask to subjective space and coregister
if true
    for i = 1:length(roiNames)
        for s = 1:length(subs)
            
            % create folder for subjective masks
            dirSubjMask = [dir.mskOut,fs,roiNames{i},'_subj',fs,'SF',num2str(subs(s),'%03d')];
            if ~exist(dirSubjMask,'dir'),mkdir(dirSubjMask),end
            cd(dirSubjMask)
            
            % copy atlas mask to subject's folder
            fileSource = [dir.mskOut,fs,'_useNow',fs,roiNames{i},'.nii'];
            copyfile(fileSource)
            
            % select inverse deformation images from T1 segmentation step
            dirStruct = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
            d = spm_select('List', dirStruct, '^iy_.*\.nii$');
            y_file = {[dirStruct fs d]}; clear d dirStruct
            job1{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            
            % select mask
            msk_file = {[dirSubjMask,fs,roiNames{i},'.nii']};
            job1{1}.spatial{1}.normalise{1}.write.subj.resample = msk_file;
            
            % defaults
            job1{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job1{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job1{1}.spatial{1}.normalise{1}.write.woptions.interp = 4; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d')])
            spm_jobman('run',job1);
            clear job1 msk_file y_file
            
            % delete unwarped file
            delete([dirSubjMask,fs,roiNames{i},'.nii'])
            
            % load sample EPI from current subject for coregistration
            dirFun = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
            d = spm_select('List', dirFun, '^uaf.*\.nii$');
            d = d(end-1,:);
            epi_file = [dirFun fs d];
            job2{1}.spm.spatial.coreg.write.ref = {epi_file};
            
            % select mask to coregister 
            job2{1}.spm.spatial.coreg.write.source = {[dirSubjMask,fs,'w',roiNames{i},'.nii']};
            
            % defaults
            job2{1}.spm.spatial.coreg.write.roptions.interp = 4;
            job2{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            job2{1}.spm.spatial.coreg.write.roptions.mask = 0;
            job2{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            
            % run job
            disp(['Coregistering sub#', num2str(subs(s),'%03d')])
            d = spm_jobman('run',job2);
            clear job2 d
            
        end
    end
end
