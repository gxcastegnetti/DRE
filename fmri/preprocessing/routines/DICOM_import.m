function DICOM_import(dir,subs,runType)
%% function DICOM_import(dirMain,subs,runType)
% ~~~
% INPUTS:
%   dirMain: root directory
%   subs: subjects
%   runtType: cell array of run types ('fun','fmap','struct','ignore')
% ~~~
% GX Castegnetti --- start ~ 05.06.18 --- last ~ 13.06.18

fs = filesep;

for s = 1:length(subs)
    
    % select run type sequence for the current subject
    runType_sub = runType{subs(s)};
    
    dirSub = [dir.dre,'data',fs,'fmri',fs,'scanner',fs,'SF',num2str(subs(s),'%03d')];
    
    if ~exist([dirSub,filesep,'fun'],'dir')
        
        % select folders with raw data (and order them properly)
        [~,dirList] = spm_select('List',[dirSub,fs,'_ima']);
        if size(dirList,1) > 9
            foo = str2num(dirList(:,end-1:end)); %#ok<*ST2NM>
        else
            foo = str2num(dirList(:,end));
        end
        [~,fooSort] = sort(foo);
        dirList = dirList(fooSort,:); clear foo fooSort
        
        % create subfolders
        funDir = [dirSub,fs,'fun'];
        structDir = [dirSub,fs,'struct'];
        fmapDir = [dirSub,fs,'fmap'];
        mkdir(funDir)
        mkdir(structDir)
        mkdir(fmapDir)
        
        % run over runs
        idxFun = 0;
        idxFmap = 0;
        for r = 1:length(runType_sub)
            
            % ignore useless runs
            if strcmp(runType_sub(r),'ignore') || strcmp(runType_sub(r),'loc')
                continue
            end
            
            % list of files
            imaDir = [dirSub,fs,'_ima',fs,dirList(r,:)];
            imaDir = imaDir(find(~isspace(imaDir))); %#ok<FNDSB>
            fileList = spm_select('List',imaDir,'ima');
            fileList = cellstr([repmat([imaDir(:,:) fs],size(fileList,1),1) fileList]);
            
            % create session-specific subfolders
            if strcmp(runType_sub(r),'fun')
                idxFun = idxFun+1;
                outDir = [funDir,fs,'S',num2str(idxFun)];
                mkdir(outDir)
            end
            if strcmp(runType_sub(r),'fmapm')
                idxFmap = idxFmap+1;
                outDir = [fmapDir,fs,'S',num2str(idxFmap),fs,'magni'];
                mkdir(outDir), cd(outDir), clear fooDir
            end
            if strcmp(runType_sub(r),'fmapp')
                outDir = [fmapDir,fs,'S',num2str(idxFmap),fs,'phase'];
                mkdir(outDir), cd(outDir), clear fooDir
            end
            if strcmp(runType_sub(r),'struct')
                outDir = structDir;
            end
            
            % prepare batch
            job{1}.spm.util.import.dicom.data = fileList;
            job{1}.spm.util.import.dicom.root = 'flat';
            job{1}.spm.util.import.dicom.outdir = {outDir};
            job{1}.spm.util.import.dicom.protfilter = '.*';
            job{1}.spm.util.import.dicom.convopts.format = 'nii';
            job{1}.spm.util.import.dicom.convopts.meta = 0;
            job{1}.spm.util.import.dicom.convopts.icedims = 0;
            
            % update user
            disp(['Importing DICOM for sub#',num2str(subs(s),'%03d'),'...']);
            
            % run batch
            spm_jobman('run',job);
            clear job
            
        end
    end
end
