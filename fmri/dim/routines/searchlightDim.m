function smm_t2 = searchlightDim(fullBrainVolumes, mask, userOptions, localOptions)

% ARGUMENTS
% fullBrainVolumes	A voxel x condition x session matrix of activity
% 				patterns.
%
% mask     		A 3d or 4d mask to perform the searchlight in.
%
% userOptions and localOptions
%
% RETURN VALUES
% smm_rs        4D array of 3D maps (x by y by z by model index) of
%               correlations between the searchlight pattern similarity
%               matrix and each of the model similarity matrices.
%
% smm_ps        4D array of 3D maps (x by y by z by model index) of p
%               values computed for each corresponding entry of smm_rs.
%
% n             an array of the same dimensions as the volume, which
%               indicates for each position how many voxels contributed
%               data to the corresponding values of the infomaps.
%               this is the number of searchlight voxels, except at the
%               fringes, where the searchlight may illuminate voxels
%               outside the input-data mask or voxel with all-zero
%               time-courses (as can arise from head-motion correction).
%
% mappingMask_actual
%               3D mask indicating locations for which valid searchlight
%               statistics have been computed.
%
% GX Castegnetti --- 2018
% Based on Niko Kriegeskorte's searchlightMapping_RDMs.m

localOptions = setIfUnset(localOptions, 'averageSessions', true);

%% here I need to take three different t_pats: one for each of the first 3 sessions
t_s1 = fullBrainVolumes(:,:,1)';
t_s2 = fullBrainVolumes(:,:,2)';
t_s3 = fullBrainVolumes(:,:,3)';
t_s4 = fullBrainVolumes(:,:,4)';
clear fullBrainVolumes;

%% Get parameters
voxSize_mm  = userOptions.voxelSize;
rad_mm      = userOptions.searchlightRadius;
monitor     = localOptions.monitor;

%% Prepare mask
inputDataMask       = logical(mask);
mappingMask_request = logical(mask);

% check to see if there's more data than mask
if (size(t_s1,2)>sum(inputDataMask(:)))
    t_s1 = t_s1(:,inputDataMask(:));
    t_s2 = t_s2(:,inputDataMask(:));
    t_s3 = t_s3(:,inputDataMask(:));
    t_s4 = t_s4(:,inputDataMask(:));
end

% find radius in voxels
volSize_vox   = size(inputDataMask);
rad_vox       = rad_mm./voxSize_mm;
minMargin_vox = floor(rad_vox);

%% create spherical multivariate searchlight
[x,y,z] = meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
sphere  = ((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(rad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
sphereSize_vox = [size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )

if monitor, figure(50); clf; showVoxObj(sphere); end % show searchlight in 3D

% compute center-relative sphere SUBindices
[sphereSUBx,sphereSUBy,sphereSUBz] = ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
sphereSUBs = [sphereSUBx,sphereSUBy,sphereSUBz];
ctrSUB     = sphereSize_vox/2+[.5 .5 .5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
ctrRelSphereSUBs = sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices

%% define masks
validInputDataMask = inputDataMask;

sumAbsY = sum(abs(t_s1),1);

validYspace_logical = (sumAbsY ~= 0) & ~isnan(sumAbsY); clear sumAbsY;
validInputDataMask(inputDataMask) = validYspace_logical; % define valid-input-data brain mask

% reduce t_pats to the valid-input-data brain mask
t_s1 = t_s1(:,validYspace_logical);
t_s2 = t_s2(:,validYspace_logical);
t_s3 = t_s3(:,validYspace_logical);
t_s4 = t_s4(:,validYspace_logical);
nVox_validInputData = size(t_s1,2);

mappingMask_request_INDs = find(mappingMask_request);
nVox_mappingMask_request = length(mappingMask_request_INDs);

if monitor
    disp([num2str(round(nVox_mappingMask_request/prod(volSize_vox)*10000)/100),'% of the cuboid volume requested to be mapped.']);
    disp([num2str(round(nVox_validInputData/prod(volSize_vox)*10000)/100),'% of the cuboid volume to be used as input data.']);
    disp([num2str(nVox_validInputData),' of ',num2str(sum(inputDataMask(:))),' declared input-data voxels included in the analysis.']);
end

volIND2YspaceIND = nan(volSize_vox);
volIND2YspaceIND(validInputDataMask) = 1:nVox_validInputData;

% n voxels contributing to infobased t at each location
n = nan(volSize_vox);

%% similarity-graph-map the volume with the searchlight
smm_t2 = nan(volSize_vox);

if monitor
    h_progressMonitor=progressMonitor(1, nVox_mappingMask_request,  'Similarity-graph-mapping...');
end

%% THE BIG LOOP! %%

for cMappingVoxI = 1:nVox_mappingMask_request
    
    if rand < 0.001
        disp(['Computing sphere around voxel ',num2str(cMappingVoxI),' of ',num2str(nVox_mappingMask_request)])
    end
    
    [x y z] = ind2sub(volSize_vox,mappingMask_request_INDs(cMappingVoxI));
    
    % compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVoxSUBs = repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1]) + ctrRelSphereSUBs;
    
    % exclude out-of-volume voxels
    outOfVolIs = (cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
        cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
        cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));
    
    cIllVoxSUBs = cIllVoxSUBs(~outOfVolIs,:);
    
    % list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVox_volINDs = sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));
    
    % restrict searchlight to voxels inside validDataBrainMask
    cIllValidVox_volINDs = cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
    cIllValidVox_YspaceINDs = volIND2YspaceIND(cIllValidVox_volINDs);
    
    % note how many voxels contributed to this locally multivariate stat
    n(x,y,z) = length(cIllValidVox_YspaceINDs);
    if n(x,y,z) < 60, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
    
    % retain only 60 random voxels int he sphere
    foo = randperm(length(cIllValidVox_YspaceINDs));
    cIllValidVox_YspaceINDs = sort(cIllValidVox_YspaceINDs(foo(1:60)));
    
    [~,~,~,~,explVar_vec_s1,~] = pca(t_s1(:,cIllValidVox_YspaceINDs)');
    [~,~,~,~,explVar_vec_s2,~] = pca(t_s2(:,cIllValidVox_YspaceINDs)');
    [~,~,~,~,explVar_vec_s3,~] = pca(t_s3(:,cIllValidVox_YspaceINDs)');
    [~,~,~,~,explVar_vec_s4,~] = pca(t_s4(:,cIllValidVox_YspaceINDs)');
    
    explVar_soFar_s1 = cumsum(explVar_vec_s1);
    idxThr = find(explVar_soFar_s1 > localOptions.explThreshold,1);
    explVar_sess(1) = idxThr;
    
    explVar_soFar_s2 = cumsum(explVar_vec_s2);
    idxThr = find(explVar_soFar_s2 > localOptions.explThreshold,1);
    explVar_sess(2) = idxThr;
    
    explVar_soFar_s3 = cumsum(explVar_vec_s3);
    idxThr = find(explVar_soFar_s3 > localOptions.explThreshold,1);
    explVar_sess(3) = idxThr;
    
    explVar_soFar_s4 = cumsum(explVar_vec_s4);
    idxThr = find(explVar_soFar_s4 > localOptions.explThreshold,1);
    explVar_sess(4) = idxThr;
    
    % take difference of correlations across context and within context
    if mean(explVar_sess) < 1, keyboard, end
    smm_t2(x,y,z) = mean(explVar_sess);
    
end

end
