function geomDiff = searchlightMapping_geom(fullBrainVolumes, mask, userOptions, localOptions)

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
t_F = fullBrainVolumes(:,1:120)';
t_B = fullBrainVolumes(:,121:end)';
clear fullBrainVolumes;

%% Get parameters
voxSize_mm  = userOptions.voxelSize;
rad_mm      = userOptions.searchlightRadius;
monitor     = localOptions.monitor;

%% Prepare mask
inputDataMask       = logical(mask);
mappingMask_request = logical(mask);

% check to see if there's more data than mask
if (size(t_F,2)>sum(inputDataMask(:)))
    t_F = t_F(:,inputDataMask(:));
    t_B = t_B(:,inputDataMask(:));
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

sumAbsY = sum(abs(t_F),1);

validYspace_logical = (sumAbsY ~= 0) & ~isnan(sumAbsY); clear sumAbsY;
validInputDataMask(inputDataMask) = validYspace_logical; % define valid-input-data brain mask

% reduce t_pats to the valid-input-data brain mask
t_F = t_F(:,validYspace_logical);
t_B = t_B(:,validYspace_logical);
nVox_validInputData = size(t_F,2);

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
geomDiff = nan([size(mask),7140]);
for cMappingVoxI = 1:nVox_mappingMask_request
    
    if mod(cMappingVoxI,1000)==0
        if monitor
            progressMonitor(cMappingVoxI, nVox_mappingMask_request, 'Searchlight mapping Mahalanobis distance...', h_progressMonitor);
        else
            fprintf('.');
        end
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
    
    if n(x,y,z) < 2, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
    
    searchlightRDM_F = corr(t_F(:,cIllValidVox_YspaceINDs)');
%     searchlightRDM_B = corr(t_B(:,cIllValidVox_YspaceINDs)');
    
    % take difference of correlations across context and within context
%     geomDiff(x,y,z,:) = abs(single(vectorizeRDM(searchlightRDM_F - searchlightRDM_B)));
    geomDiff(x,y,z,:) = single(vectorizeRDM(searchlightRDM_F));
end

%% END OF THE BIG LOOP! %%

if monitor
    fprintf('\n');
    close(h_progressMonitor);
end

mappingMask_actual = mappingMask_request;
mappingMask_actual(isnan(sum(smm_t2,4))) = 0;

end