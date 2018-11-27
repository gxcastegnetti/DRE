function userOptions = dre_rsa_userOptions(dir,subs)

% subjects
for s = 1:length(subs)
    userOptions.subjectNames{s} = ['SF',num2str(subs(s),'%03d')];
end

% default colour label for RDMs corresponding to RoI masks (as opposed to models)
userOptions.RoIColor = [0 0 1];
userOptions.ModelColor = [0 1 0];

%% First-order analysis

% read object names
objVersion  = 7; % set which column to read according to the object set version
objs        = readtable([dir.behDat,filesep,'Objects.csv']);
foo         = logical(table2array(objs(:,objVersion))); 
objsName    = table2cell(objs(foo,2)); clear foo objs

% attach F and B conditions
for h = 1:length(objsName)
    objsF{h} = ['F-',objsName{h}];
    objsB{h} = ['B-',objsName{h}];
end
condLabels = [objsF';objsB'];

% text lables which may be attached to the conditions for MDS plots
userOptions.conditionLabels = condLabels;

% colours for the conditions
userOptions.conditionColours = kron([1 0 0; 0 0 1], ones(120,1));

% groups to plot convex hulls around
userOptions.convexHulls = [ones(1,16) 2*ones(1,16) 3*ones(1,16) 4*ones(1,16)];

% metric to use when calculating first-order RDMs
userOptions.distance = 'mahalanobis';

%% Second-order analysis

% metric for second-level comparison
userOptions.distanceMeasure = 'Spearman';

% how many permutations should be used to test the significance of the fits (10000 highly recommended)
userOptions.significanceTestPermutations = 10000;

% should RDMs entries be rank transformed into [0,1]?
userOptions.rankTransform = true;

% RDM colourscheme
userOptions.colourScheme = jet(64);

% should rubber bands be shown on the MDS plot?
userOptions.rubberbands = false;

% what criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

% how should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = false;
userOptions.saveFiguresFig = false;
userOptions.saveFiguresPS = false;
userOptions.saveFiguresEps = false;

% bootstrap options
userOptions.nResamplings = 1000;
userOptions.resampleSubjects = true;
userOptions.resampleConditions = true;

% ceiling estimation
userOptions.RDMname = 'referenceRDM';
userOptions.plottingStyle = 2;

end
