function [varargout] = fMRIDataPreparation(betaCorrespondence, userOptions)
%
% fMRIDataPreparation is a function designed to take fMRI data and load it into
% the correct format for the rest of the toolbox to use.
%
% [fullBrainVols =] fMRIDataPreparation(betaCorrespondence, userOptions)
%
%       betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image.
%               Alternatively, this can be the string 'SPM', in which case the
%               SPM metadata will be used to infer this information, provided
%               that userOptions.conditionLabels is set, and the condition
%               labels are the same as those used in SPM.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.betaPath
%                        A string which contains the absolute path to the
%                        location of the beta images. It can contain the
%                        following wildcards which would be replaced as
%                        indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                                [[betaIdentifier]]
%                                        To be replaced by filenames as provided
%                                        by betaCorrespondence.
%                userOptions.conditionLabels
%                        A cell array containing the names of the conditions in
%                        this experiment. This must be set if the 'SPM' option
%                        is used.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_ImageData.mat
%                        Contains the raw beta images in a struct such that
%                        fullBrainVols.(subject) is a [nVoxels nConditions
%                        nSessions] matrix.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_fMRIDataPreparation_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
%  Cai Wingfield 11-2009 -- 1-2010, 6-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll return to the pwd when the function has finished

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('fMRIDataPreparation:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRIDataPreparation:NoRootPath', 'rootPath must be set. See help'); end%if
if ~isfield(userOptions, 'betaPath'), error('fMRIDataPreparation:NoBetaPath', 'betaPath must be set. See help'); end%if
if ~isfield(userOptions, 'subjectNames'), error('fMRIDataPreparation:NoSubjectNames', 'subjectNames must be set. See help'); end%if
if (~isfield(userOptions, 'conditionLabels') && ischar(betaCorrespondence) && strcmpi(betaCorrespondence, 'SPM')), error('fMRIDataPreparation:NoConditionLabels', 'conditionLables must be set if the data is being extracted from SPM.'); end%if

%% Get Data
nSubjects = numel(userOptions.subjectNames);

fprintf('Gathering scans.\n');

for subject = 1:nSubjects
    
    betas = spm.getDataFromSPM(userOptions,subject);
    nConditions = size(betas, 2);
    nSessions = size(betas, 1);
    
    % Figure out the subject's name
    thisSubject = userOptions.subjectNames{subject};
    
    fprintf(['Reading beta volumes for subject number ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject]);
    
    for session = 1:nSessions
        for condition = 1:nConditions
            
            readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(session,condition).identifier, '[[subjectName]]', thisSubject);
            if strcmp(betaCorrespondence,'SPM')
                brainMatrix = spm_read_vols(spm_vol(readPath));
            else
                load(readPath);
                brainMatrix = betaImage;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GXC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14.08.18
            brainVector = reshape(brainMatrix, 1, []);
            brainVector_full = brainVector;
            brainVector = brainVector(~isnan(brainVector));
            
            subjectMatrix(:, condition, session) = brainVector; % (voxel, condition, session)
            subjectMatrix_full(:, condition, session) = brainVector_full;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GXC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14.08.18
            
            clear brainMatrix brainVector brainVector_full;
            
            fprintf('.');
            
        end
        
    end
    
    % For each subject, record the vectorised brain scan in a subject-name-indexed structure
    fullBrainVols.(thisSubject) = subjectMatrix; clear subjectMatrix;
    fullBrainVols_full.(thisSubject) = subjectMatrix_full; clear subjectMatrix_full;
    
    fprintf('\b:\n');
    
end

varargout{1} = fullBrainVols;
varargout{2} = fullBrainVols_full;

cd(returnHere); % Go back (probably will never have left)

end%
