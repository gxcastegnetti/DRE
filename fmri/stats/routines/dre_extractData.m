function bData = dre_extractData(dir,subs,taskOrd,excNan)
%% function DRE_extractData(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dirSub: directory of subject's data
%   subj: subject number
%   taskOrd: FBFB or BFBF
% OUTPUTS:
%   bData: structure with all trial onsets and behavioural parameters
% ~~~
% GX Castegnetti --- start ~ 13.06.18 --- last ~ 22.06.18

fs = filesep;
n_sess = 4;

dirData = [dir.dre,fs,'data'];

%% read objects
objs        = readtable([dir.beh,fs,'Objects.csv']);
objsName    = table2cell(objs(:,2)); clear foo objs

%% loop over subjects and sessions
for s = 1:length(subs)
    for r = 1:n_sess
        
        %% extract subject-specific task order
        if taskOrd(s) == 1
            day2Order = {'F','B','F','B'};
            day1Order = {'1','2','1','2'};
        elseif taskOrd(s) == 2
            day2Order = {'B','F','B','F'};
            day1Order = {'1','2','1','2'};
        end
        
        %% subject directories
        dirPsyO = [dirData,fs,'fmri',fs,'psychOut',fs,'SF',num2str(subs(s),'%03d')];
        dirBeha = [dirData,fs,'behaviour',fs,'SF',num2str(subs(s),'%03d')];
        
        %% trial onsets
        
        % load behavioural matrix from day 2
        Mday2 = csvread([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'.csv']);
        
        % load time flags
        t_foo = csvread([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'_time.csv']);
        t0 = t_foo(3); clear t_foo
        
        % find imagination and choice onsets
        onsIma = Mday2(Mday2(:,3)==0,2) - t0; % imagination
        onsCho = Mday2(Mday2(:,3)==1,2) - t0; % choice
        
        %% object IDs
        
        % extract IDs from day 2
        idIma_day2 = Mday2(Mday2(:,3)==0,4);   % imagination
        idCho_day2 = Mday2(Mday2(:,3)==1,4:5); % choice
        chosenChoice_day2 = Mday2(Mday2(:,3)==1,7);
        latencChoice_day2 = Mday2(Mday2(:,3)==1,8);
        
        % object names
        objSes = objsName(idIma_day2);
        
        % load behaviour matrix from day 1
        Mday1 = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B',day1Order{r},'_DRE.csv']);
        
        %% find value that on day 1 was assigned to objects presented on day 2
        
        %% imagination
        idxIma_day2 = NaN(length(idIma_day2),1);
        for i = 1:length(idIma_day2)
            idxIma_day2(i) = find(Mday1(:,2) == idIma_day2(i));
        end
        objVal = Mday1(idxIma_day2,3); % value
        objCon = Mday1(idxIma_day2,4); % confidence
        
        %% choice
        idxCho_day2_L = NaN(length(idCho_day2),1);
        idxCho_day2_R = NaN(length(idCho_day2),1);
        choVal = NaN(length(idCho_day2),1);
        for i = 1:length(idCho_day2)
            idxCho_day2_L(i) = find(Mday1(:,2) == idCho_day2(i,1));
            idxCho_day2_R(i) = find(Mday1(:,2) == idCho_day2(i,2));
            if chosenChoice_day2(i) == -1
                choVal(i) = Mday1(idxCho_day2_L(i),3);
            elseif chosenChoice_day2(i) == 1
                choVal(i) = Mday1(idxCho_day2_R(i),3);
            else
                choVal(i) = NaN;
            end
        end
        difVal = abs(Mday1(idxCho_day2_L,3) - Mday1(idxCho_day2_R,3));
        
        % movement onset and side
        movCho = Mday2(Mday2(:,3)==1,8);
        sidCho = Mday2(Mday2(:,3)==1,7);
        
        % familiarity and value
        Mday1F = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']);
        
        % imagination
        for i = 1:length(idIma_day2)
            idxFam_day2(i) = find(Mday1F(:,2)==idIma_day2(i)); %#ok<*AGROW>
        end
        objFam = Mday1F(idxFam_day2,3); % familiarity
        objMon = Mday1F(idxFam_day2,4); % monetary value
        
        % remove NaNs
        if excNan == 1
            objKeep = ~isnan(objVal) & ~isnan(objCon) & ~isnan(objFam) & ~isnan(objMon);
            difKeep = ~isnan(difVal) & ~isnan(choVal) & ~isnan(movCho) & ~isnan(sidCho);
            onsIma = onsIma(objKeep);
            objVal = objVal(objKeep);
            objCon = objCon(objKeep);
            objFam = objFam(objKeep);
            objMon = objMon(objKeep);
            objSes = objSes(objKeep);
            onsCho = onsCho(difKeep);
            difVal = difVal(difKeep);
            choVal = choVal(difKeep);
            movCho = movCho(difKeep);
            sidCho = sidCho(difKeep);
        end
        
        if strcmp(day2Order{r},'F')
            
            bData(subs(s)).sessType{r} = 'fire';
            
            % fire
            bData(subs(s)).imagination(r).fire.names = objSes;
            bData(subs(s)).imagination(r).fire.valIdx = idxIma_day2;
            bData(subs(s)).imagination(r).fire.famIdx = idxFam_day2;
            bData(subs(s)).imagination(r).fire.onset = onsIma;
            bData(subs(s)).imagination(r).fire.value = objVal;
            bData(subs(s)).imagination(r).fire.confidence = objCon;
            bData(subs(s)).imagination(r).fire.familiarity = objFam;
            bData(subs(s)).imagination(r).fire.econValue = objMon;
            bData(subs(s)).choice(r).fire.onset = onsCho;
            bData(subs(s)).choice(r).fire.valueDiff = difVal;
            bData(subs(s)).choice(r).fire.valueChosen = choVal;
            bData(subs(s)).choice(r).fire.movOnset = movCho;
            bData(subs(s)).choice(r).fire.movSide = sidCho;
            
            % boat
            bData(subs(s)).imagination(r).boat.onset = [];
            bData(subs(s)).imagination(r).boat.value = [];
            bData(subs(s)).imagination(r).boat.confidence = [];
            bData(subs(s)).imagination(r).boat.familiarity = [];
            bData(subs(s)).imagination(r).boat.econValue = [];
            bData(subs(s)).choice(r).boat.onset = [];
            bData(subs(s)).choice(r).boat.valueDiff = [];
            bData(subs(s)).choice(r).boat.valueChosen = [];
            bData(subs(s)).choice(r).boat.movOnset = [];
            bData(subs(s)).choice(r).boat.movSide = [];
            
        elseif strcmp(day2Order{r},'B')
            
            bData(subs(s)).sessType{r} = 'boat';
            
            % fire
            bData(subs(s)).imagination(r).fire.onset = [];
            bData(subs(s)).imagination(r).fire.value = [];
            bData(subs(s)).imagination(r).fire.confidence = [];
            bData(subs(s)).imagination(r).fire.familiarity = [];
            bData(subs(s)).imagination(r).fire.econValue = [];
            bData(subs(s)).choice(r).fire.onset = [];
            bData(subs(s)).choice(r).fire.valueDiff = [];
            bData(subs(s)).choice(r).fire.valueChosen = [];
            bData(subs(s)).choice(r).fire.movOnset = [];
            bData(subs(s)).choice(r).fire.movSide = [];
            
            % boat
            bData(subs(s)).imagination(r).boat.names = objSes;
            bData(subs(s)).imagination(r).boat.onset = onsIma;
            bData(subs(s)).imagination(r).boat.valIdx = idxIma_day2;
            bData(subs(s)).imagination(r).boat.famIdx = idxFam_day2;
            bData(subs(s)).imagination(r).boat.value = objVal;
            bData(subs(s)).imagination(r).boat.confidence = objCon;
            bData(subs(s)).imagination(r).boat.familiarity = objFam;
            bData(subs(s)).imagination(r).boat.econValue = objMon;
            bData(subs(s)).choice(r).boat.onset = onsCho;
            bData(subs(s)).choice(r).boat.valueDiff = difVal;
            bData(subs(s)).choice(r).boat.valueChosen = choVal;
            bData(subs(s)).choice(r).boat.movOnset = movCho;
            bData(subs(s)).choice(r).boat.movSide = sidCho;
        end
        
    end
end