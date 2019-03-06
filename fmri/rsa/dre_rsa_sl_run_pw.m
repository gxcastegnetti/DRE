%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_ima';
betaid       = 'rsa_pulse_ima';

%% directories
fs         = filesep;
dir.rsaCod = pwd;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.betaid = betaid;

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 7 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39:43 47:49];
taskOrd = [ones(1,10),2*ones(1,10),1,2,ones(1,5),2*ones(1,4)];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% 1st level
roiName = 'none';
if false
    nameBeta = ['level1',fs,betaid,fs,roiName];
    bData = dre_extractData(dir,subs,taskOrd,0);
    timing.iOns = 0;
    timing.iDur = 0;
    dre_level1_rsa(dir,nameBeta,subs,bData,timing,roiName);
end

%% extract models of value, confidence, familiarity, price
RDMs = dre_extractRDMs(dir,subs,taskOrd);

%% run searchlight for imagination
% create matrix for object ID
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];

%% materials and weights
dirMaterials = [dir.behDat,filesep,'materials'];
est_MZ = csvread([dirMaterials,filesep,'ObjectsMaterials_MZ.csv'],1,2);
est_LU = csvread([dirMaterials,filesep,'ObjectsMaterials_LU.csv'],1,2);
est_ND = csvread([dirMaterials,filesep,'ObjectsMaterials_ND.csv'],1,2);
est_LU = est_LU(est_LU(:,5) ~= 0,:);
est_ND = est_ND(est_ND(:,5) ~= 0,:);
est_mean = (est_MZ + est_LU + est_ND)/3;

weight = [est_mean(:,5);est_mean(:,5)];
wood   = [est_mean(:,1);est_mean(:,1)];
metal  = [est_mean(:,2);est_mean(:,2)];
plastic = [est_mean(:,3);est_mean(:,3)];
fabric = [est_mean(:,4);est_mean(:,4)];
for i = 1:numel(weight)
    for j = 1:numel(weight)
        RDM_weight(i,j) = abs(weight(i) - weight(j));
        RDM_material(i,j) = sqrt((wood(i) - wood(j))^2 + (metal(i) - metal(j))^2 + ...
            (plastic(i) - plastic(j))^2 + (fabric(i) - fabric(j))^2);
    end
end, clear weight wood metal plastic fabric

%% compute correlation between value and weight/material and between goals
RDM_weight_vec = vectorizeRDM(RDM_weight);

foo = eye(120);
foo(foo == 1) = nan;

fooMat = [zeros(120),foo;foo,zeros(120)];

RDM_material_vec = vectorizeRDM(RDM_material + fooMat);
for s = 1:length(subs)
    RDM_value_vec = vectorizeRDM(RDMs{s}.val);
    [r_valVSmat(s),p_valVSmat(s)] = corr(RDM_value_vec',RDM_material_vec','type','spearman','rows','pairwise');
    [r_valVSwei(s),p_valVSwei(s)] = corr(RDM_value_vec',RDM_weight_vec','type','spearman','rows','pairwise');
    
    % goals
    [r_goals(s),p_goals(s)] = corr(vectorizeRDM(RDMs{s}.val_F)',vectorizeRDM(RDMs{s}.val_B)','type','spearman','rows','pairwise');
end


%% searchlight
for s = 1:length(subs)
    
    RDMs{s}.oid = 1 - mat_ID;
    
    % update user
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    
    % prepare mask
    fileMask = [dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
    
    %% run searchlight
    model = RDM_material_vec;
    rs = searchlight_pw(dir,subs(s),analysisName,fileMask,model); %#ok<*ASGLU>
    save([dir.out,fs,analysisName,fs,'sl_valPcon_SF',num2str(subs(s),'%03d')],'rs','model')
    clear model rs binaryMask
end
