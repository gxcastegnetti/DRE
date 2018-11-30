function rs = searchlight_pw(dir,subNum,analysisName,fileMask,model)

% this is an example script for doing RSA searchlight from the outputs of
% SPM's first-level analysis.
% you will need to change the paths, etc.
% nili.hamed@gmail.com

fs = filesep;

%% obtain the SL definitions
% volumetric SL
Vmask = spm_vol(fileMask);
Vmask.data = spm_read_vols(Vmask);
L = defineSearchlight({Vmask},Vmask,'sphere',9);
fname = [dir.out,fs,analysisName,fs,'sl_SF',num2str(subNum,'%03d')];
save(fname,'-struct','L');

%% run the SL
% load SL definition
L = load(fname);
subjDir = [dir.rsaOut,fs,'level1',fs,dir.betaid,fs,'none',fs,'SF',num2str(subNum,'%03d')];
cd(subjDir);

load('SPM.mat')

% Define output images
outFiles={};
for r = 1:1 % number of runs
    for k = 1:28680 % number of dissims
        outFiles{end+1} = sprintf('dists_run%d.nii,%d',r,k); % these are the name of the files that will be written in the subject's folder (where the SPM.mat of your subject is)
    end
end

% this part of the code extract the indices of the conditions of interest
% from the design matrix.
% for r = 1:4
%     vec = zeros(1,size(SPM.xX.xKXs.X,2));
%     u = vec;
%     idx = SPM.Sess(r).col(1:60);
%     vec(idx) = 1:60;
%     condVecs{r} = vec;
% end
condVecs{1} = [true(1,60),false(1,6),true(1,60),false(1,6),true(1,60),false(1,6),true(1,60),false(1,10)];
rs = runSearchlight(L,SPM.xY.VY,outFiles,model,subNum,@rsaFunc,'optionalParams',{SPM,condVecs});

function out = rsaFunc(Y,SPM,condVecs)
B = noiseNormalizeBeta(Y,SPM);      % Get prewhitened beta weights
for r = 1:1
    thisB = B(condVecs{r},:);
    RDM = pdist(real(thisB),'correlation');
    out = RDM(:);     % Arrange the outputs in a vector, as they are written to files
end
