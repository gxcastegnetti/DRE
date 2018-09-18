function spatial_preprocess(dir,subs,fmapFiles)
% function spatialPreprocessing(dirSPM,dirSub,subs,fmapFiles)
% ~~~
% This script:
% a) corrects slice timing
% b) realigns & unwarps using preprocessed fieldmaps
% c) coregisters the structural to mean-epi generated in (b)
% d) segments the coregistered structural
% e) normalizes all functional runs using forward deformations generated in (d)...
% f) ...or normalises structural to MNI template and applies warp to functional runs
% g) smooths functional images
% ~~~
% GX Castegnetti --- start ~ 05.06.18 --- last ~ 31.07.18 (adapted from SM Fleming's code)

fs = filesep;

%% parameters
slicetiming = 1;
realign     = 1;
coregister  = 1;
smoothNat   = 1;
segment     = 1;
normalise   = 1;
smooth      = 1;
FWHM        = 8;
resolEPI    = [2 2 2]; %#ok<*NASGU>
nslices     = 48;
TR          = 3.36;
sliceorder  = 1:nslices;

%% loop
for s = 1:length(subs)
    
    n_sess = length(fmapFiles{subs(s)});
    
    % subject directory
    dirSub = [dir.dre,'data/fmri/scanner/SF',num2str(subs(s),'%03d')];
    
    %% Slice-timing
    if slicetiming
        disp(['Slice timing job specification for sub#', num2str(subs(s),'%03d'),'...']);
        
        % loop over sessions
        for r = 1:n_sess
            % functional directory
            dirFun = [dirSub,fs,'fun',fs,'S',num2str(r)];
            
            % select epi files and assign to job
            d = spm_select('List', dirFun, '^f.*\.nii$');
            files = cellstr([repmat([dirFun fs],size(d,1),1) d]);
            job{1}.spm.temporal.st.scans(r) = {files}; clear d files dirSess
        end
        job{1}.spm.temporal.st.nslices = nslices;
        job{1}.spm.temporal.st.tr = TR;
        job{1}.spm.temporal.st.ta = TR - (TR/nslices);
        job{1}.spm.temporal.st.so = sliceorder;
        job{1}.spm.temporal.st.refslice = round(nslices/2);
        job{1}.spm.temporal.st.prefix = 'a';
        
        % save and run job
        disp(['Correcting slice timing for sub#',num2str(subs(s),'%03d'),'...']);
        d = spm_jobman('run',job);
        clear job
    end
    
    
    %% Realignment and unwarping
    if realign
        disp(['Realign job specification for sub#',num2str(subs(s),'%03d'),'...']);
        cd(dirSub)
        % loop over sessions
        for r = 1:n_sess
            % select directory
            dirFun = [dirSub,'/fun/S',num2str(r)];
            
            % select scans and assign to job
            if slicetiming
                d = spm_select('List', dirFun, '^af.*\.nii$');
            else
                d = spm_select('List', dirFun, '^f.*\.nii$'); %#ok<*UNRCH>
            end
            files  = cellstr([repmat([dirFun fs],size(d,1),1) d]);
            job{1}.spatial{1}.realignunwarp.data(r).scans = files; clear f files
            
            % select fieldmap aligned to first epi in this session, add to structure
            fmfilename = fmapFiles{subs(s)}{r};
            job{1}.spatial{1}.realignunwarp.data(r).pmscan{1} = fmfilename;   % Note address the cell first with {}, then the structure with (), then put in cell...
        end
        
        % specify estimation options
        job{1}.spatial{1}.realignunwarp.eoptions.quality = 0.9;      % default 0.9. Steve 1
        job{1}.spatial{1}.realignunwarp.eoptions.sep     = 4.0;
        job{1}.spatial{1}.realignunwarp.eoptions.fwhm    = 5.0;
        job{1}.spatial{1}.realignunwarp.eoptions.rtm     = 0;
        job{1}.spatial{1}.realignunwarp.eoptions.einterp = 2;        % default 2, Steve 7
        job{1}.spatial{1}.realignunwarp.eoptions.ewrap   = [0 0 0];
        job{1}.spatial{1}.realignunwarp.eoptions.weight  = {};
        
        job{1}.spatial{1}.realignunwarp.uweoptions.basfcn   = [12 12];
        job{1}.spatial{1}.realignunwarp.uweoptions.regorder = 1.0;
        job{1}.spatial{1}.realignunwarp.uweoptions.lambda   = 1e+005;
        job{1}.spatial{1}.realignunwarp.uweoptions.jm       = 0;
        job{1}.spatial{1}.realignunwarp.uweoptions.fot      = [4 5];
        job{1}.spatial{1}.realignunwarp.uweoptions.sot      = [];
        job{1}.spatial{1}.realignunwarp.uweoptions.uwfwhm   = 4;
        job{1}.spatial{1}.realignunwarp.uweoptions.rem      = 1;
        job{1}.spatial{1}.realignunwarp.uweoptions.noi      = 5;
        job{1}.spatial{1}.realignunwarp.uweoptions.expround = 'Average';
        
        % specify reslice options
        job{1}.spatial{1}.realignunwarp.uwroptions.uwwhich = [2 1];
        job{1}.spatial{1}.realignunwarp.uwroptions.rinterp = 4;      % default 4, Steve 7
        job{1}.spatial{1}.realignunwarp.uwroptions.wrap    = [0 0 0];
        job{1}.spatial{1}.realignunwarp.uwroptions.mask    = 1.0;
        job{1}.spatial{1}.realignunwarp.uwroptions.prefix  = 'u';
        
        % save and run job
        disp(['Realigning sub#', num2str(subs(s),'%03d'),'...']);
        d = spm_jobman('run',job);
        clear job
    end
    
    
    %% Coregister
    if coregister
        cd(dirSub)
        % coregister the structural to the mean epi from the first run
        disp(['Coregister job specification for sub#', num2str(subs(s),'%03d'),'...']);
        
        % select source directory (select mean epi as reference)
        dirSource = [dirSub,fs,'fun',fs,'S1'];
        d = spm_select('List', dirSource, '^meanua.*\.nii$');
        
        files = [dirSource fs d];
        job{1}.spm.spatial.coreg.estimate.ref = {files};
        
        % specify structural ref image
        dirStruct = [dirSub,fs,'struct'];
        d = spm_select('List', dirStruct, '^s.*\.nii$');
        fileStruct = [dirStruct fs d];
        job{1}.spm.spatial.coreg.estimate.source            = {fileStruct}; clear fileStruct
        job{1}.spm.spatial.coreg.estimate.other             = {''};
        job{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        job{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
        job{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.0200 0.0200 0.0200 0.0010 0.0010 0.0010 ...
            0.0100 0.0100 0.0100 0.0010 0.0010 0.0010];
        job{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
        job{1}.spm.spatial.coreg.estwrite.roptions.interp   = 4;
        job{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
        job{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
        job{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'c';
        
        % run job
        disp(['Coregistering sub#',num2str(subs(s),'%03d'),'...']);
        d = spm_jobman('run',job);
        clear job d files
    end
    
    %% light smoothing on native space (for RSA)
    if smoothNat
        disp(['Native smoothing job specification for sub#', num2str(subs(s),'%03d')]);
        
        % get normalized scans
        conCat_files = [];
        for r = 1:n_sess
            dirFun = [dirSub,fs,'fun',fs,'S',num2str(r)];
            d = spm_select('List', dirFun, '^uaf.*\.nii$');
            files  = cellstr([repmat([dirFun fs],size(d,1),1) d]);
            conCat_files = [conCat_files; files]; clear d files % concatenate all files over runs
        end
        
        job{1}.spatial{1}.smooth.data                          = conCat_files;
        job{1}.spatial{1}.smooth.fwhm                          = [3 3 3]; % <--- light smoothing
        job{1}.spatial{1}.smooth.dtype                         = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix               = 's';
        
        % save and run job
        disp(['Smoothing sub#',num2str(subs(s),'%03d'),'...'])
        d = spm_jobman('run',job);
        clear job
    end
    
    %% Unified segmentation (new segment)
    if segment
        disp(['Segmentation job specification for sub#', num2str(subs(s),'%03d'),'...']);
        
        dirStruct = [dirSub,fs,'struct'];
        d = spm_select('List', dirStruct, '^s.*\.nii$'); % notice coreg only changed header .mat file, so same address
        fileStruct = {[dirStruct fs d]};
        
        job{1}.spatial{1}.preproc.channel.vols = fileStruct;
        job{1}.spatial{1}.preproc.channel.biasreg = 1e-03;
        job{1}.spatial{1}.preproc.channel.biasfwhm = 60;
        job{1}.spatial{1}.preproc.channel.write = [1 0];   % save bias-corrected image
        
        % Tissue priors
        job{1}.spatial{1}.preproc.tissue(1).tpm{1} = [dir.spm fs 'tpm/TPM.nii,1'];
        job{1}.spatial{1}.preproc.tissue(1).ngaus = 1;
        job{1}.spatial{1}.preproc.tissue(1).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(1).warped = [0 0];
        job{1}.spatial{1}.preproc.tissue(2).tpm{1} = [dir.spm fs 'tpm/TPM.nii,2'];
        job{1}.spatial{1}.preproc.tissue(2).ngaus = 1;
        job{1}.spatial{1}.preproc.tissue(2).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(2).warped = [0 0];
        job{1}.spatial{1}.preproc.tissue(3).tpm{1} = [dir.spm fs 'tpm/TPM.nii,3'];
        job{1}.spatial{1}.preproc.tissue(3).ngaus = 2;
        job{1}.spatial{1}.preproc.tissue(3).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(3).warped = [0 0];
        job{1}.spatial{1}.preproc.tissue(4).tpm{1} = [dir.spm fs 'tpm/TPM.nii,4'];
        job{1}.spatial{1}.preproc.tissue(4).ngaus = 3;
        job{1}.spatial{1}.preproc.tissue(4).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(4).warped = [0 0];
        job{1}.spatial{1}.preproc.tissue(5).tpm{1} = [dir.spm fs 'tpm/TPM.nii,5'];
        job{1}.spatial{1}.preproc.tissue(5).ngaus = 4;
        job{1}.spatial{1}.preproc.tissue(5).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(5).warped = [0 0];
        job{1}.spatial{1}.preproc.tissue(6).tpm{1} = [dir.spm fs 'tpm/TPM.nii,6'];
        job{1}.spatial{1}.preproc.tissue(6).ngaus = 2;
        job{1}.spatial{1}.preproc.tissue(6).native = [1 0];
        job{1}.spatial{1}.preproc.tissue(6).warped = [0 0];
        
        % Warp options
        job{1}.spatial{1}.preproc.warp.mrf = 1;
        job{1}.spatial{1}.preproc.warp.cleanup = 1;
        job{1}.spatial{1}.preproc.warp.reg = [0 1e-03 0.5 0.05 0.2];
        job{1}.spatial{1}.preproc.warp.affreg = 'mni';
        job{1}.spatial{1}.preproc.warp.fwhm = 0;
        job{1}.spatial{1}.preproc.warp.samp = 3;
        job{1}.spatial{1}.preproc.warp.write = [1 1];  % write out inverse + forward deformations
        
        % save and run job
        disp(['Segmenting sub#', num2str(subs(s),'%03d'),'...']);
        d = spm_jobman('run',job);
        clear job
    end
    
    %% Normalise bias-correct structural and slicetime corrected EPIs
    if normalise
        disp(['Normalisation job specification for sub#', num2str(subs(s),'%03d'),'...']);
                
        % Loop over sessions for epi's
        conCat_files = [];
        for r = 1:n_sess
            dirFun = [dirSub,fs,'fun',fs,'S',num2str(r)];
            
            % select unwarped files
            d = spm_select('List', dirFun, '^uaf.*\.nii$');
            files  = cellstr([repmat([dirFun fs],size(d,1),1) d]);
            conCat_files = [conCat_files; files]; clear d files %#ok<*AGROW> % concatenate all files over runs
        end
        
        dirStruct = [dirSub,fs,'struct'];
        d = spm_select('List', dirStruct, '^s.*\.nii$'); % notice coreg only changed header .mat file, so same address
        fileStruct = {[dirStruct fs d]};
        conCat_files = [conCat_files; fileStruct];
        
        % select forward deformation images from T1 segmentation step
        dirStruct = [dirSub,fs,'struct'];
        d = spm_select('List', dirStruct, '^y_.*\.nii$');
        files  = {[dirStruct fs d]};
        job{1}.spatial{1}.normalise{1}.write.subj.def = files; clear d files
        
        job{1}.spatial{1}.normalise{1}.write.subj.resample     = conCat_files;
        job{1}.spatial{1}.normalise{1}.write.woptions.bb       = [-78 -112 -70; 78 76 85]; % bounding box of volume
        job{1}.spatial{1}.normalise{1}.write.woptions.vox      = resolEPI;  % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
        job{1}.spatial{1}.normalise{1}.write.woptions.interp   = 4;        % changed default to 7th degree B-spline
        
        % save and run job
        disp(['Normalising sub#', num2str(subs(s),'%03d')])
        d = spm_jobman('run',job);
        clear job
    end
    
    %% Smoothing
    if smooth
        disp(['Smoothing job specification for sub#', num2str(subs(s),'%03d')]);
        
        % get normalized scans
        conCat_files = [];
        for r = 1:n_sess
            dirFun = [dirSub,fs,'fun',fs,'S',num2str(r)];
            d = spm_select('List', dirFun, '^wuaf.*\.nii$');
            files  = cellstr([repmat([dirFun fs],size(d,1),1) d]);
            conCat_files = [conCat_files; files]; clear d files % concatenate all files over runs
        end
        
        job{1}.spatial{1}.smooth.data                          = conCat_files;
        job{1}.spatial{1}.smooth.fwhm                          = [FWHM FWHM FWHM];
        job{1}.spatial{1}.smooth.dtype                         = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix               = 's';
        
        % save and run job
        disp(['Smoothing sub#',num2str(subs(s),'%03d'),'...'])
        d = spm_jobman('run',job);
        clear job
    end
    
end
