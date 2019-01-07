function bData = dre_extractData(dir,subs,taskOrd,excNan)
%% function dre_extractData(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dir: struct of directories
%   subj: subject numbers
%   taskOrd: FBFB or BFBF
%   excNaN: whether (1) or not (0) to exclude NaNs from behavioural
%           measures
% OUTPUTS:
%   bData: struct with trial onsets and behavioural measures
% ~~~
% GX Castegnetti --- 2018

fs = filesep;
n_sess = 4;

dirData = [dir.dre,fs,'data'];

%% read objects
objs        = readtable([dir.behDat,fs,'Objects.csv']);
objsNameAll = table2cell(objs(:,2)); clear foo objs

%% loop over subjects and sessions
for s = 1:length(subs)
    for r = 1:n_sess
        
        %% extract subject-specific task order
        if taskOrd(s) == 1
            day2Order = {'F','B','F','B'};
            day1Order = {'1','2','1','2'};
            day1OrderWrong = {'2','1','2','1'};
        elseif taskOrd(s) == 2
            day2Order = {'B','F','B','F'};
            day1Order = {'2','1','2','1'};
            day1OrderWrong = {'1','2','1','2'};
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
        
        % object names
        objSes = objsNameAll(idIma_day2);
        
        % load behaviour matrix from day 1
        Mday1 = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B',day1Order{r},'_DRE.csv']);
        Mday1Wrong = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B',day1OrderWrong{r},'_DRE.csv']);
        
        %% find value that on day 1 was assigned to objects presented on day 2
        
        %% imagination
        idxIma_day2 = NaN(length(idIma_day2),1);
        for i = 1:length(idIma_day2)
            idxIma_day2(i) = find(Mday1(:,2) == idIma_day2(i));
        end
        objVal = Mday1(idxIma_day2,3); % value
        objValWrong = Mday1Wrong(idxIma_day2,3); % value
        objCon = Mday1(idxIma_day2,4); % confidence
        
        Mday1F = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']);
        
        % familiarity and price
        for i = 1:length(idIma_day2)
            idxFam_day2(i) = find(Mday1F(:,2) == idIma_day2(i)); %#ok<*AGROW>
        end
        objFam = Mday1F(idxFam_day2,3); % familiarity
        objPri = Mday1F(idxFam_day2,4); % monetary value
        
        %% choice
        idxCho_day2_L = NaN(length(idCho_day2),1);
        idxCho_day2_R = NaN(length(idCho_day2),1);
        valCho = NaN(length(idCho_day2),1);
        valUnc = NaN(length(idCho_day2),1);
        for i = 1:length(idCho_day2)
            idxCho_day2_L(i) = find(Mday1(:,2) == idCho_day2(i,1));
            idxCho_day2_R(i) = find(Mday1(:,2) == idCho_day2(i,2));
            if chosenChoice_day2(i) == -1
                valCho(i) = Mday1(idxCho_day2_L(i),3);
                valUnc(i) = Mday1(idxCho_day2_R(i),3);
            elseif chosenChoice_day2(i) == 1
                valCho(i) = Mday1(idxCho_day2_R(i),3);
                valUnc(i) = Mday1(idxCho_day2_L(i),3);
            else
                valCho(i) = NaN;
                valUnc(i) = NaN;
            end
        end
        choObjPresented = [idxCho_day2_L idxCho_day2_R];
        difVal = abs(valCho - valUnc); % value difference
        chMunc = valCho - valUnc; % value chosen - value unchosen
        
        %% choice wrong
        idxCho_day2_L = NaN(length(idCho_day2),1);
        idxCho_day2_R = NaN(length(idCho_day2),1);
        valChoWrong = NaN(length(idCho_day2),1);
        valUncWrong = NaN(length(idCho_day2),1);
        for i = 1:length(idCho_day2)
            idxCho_day2_L(i) = find(Mday1Wrong(:,2) == idCho_day2(i,1));
            idxCho_day2_R(i) = find(Mday1Wrong(:,2) == idCho_day2(i,2));
            if chosenChoice_day2(i) == -1
                valChoWrong(i) = Mday1(idxCho_day2_L(i),3);
                valUncWrong(i) = Mday1(idxCho_day2_R(i),3);
            elseif chosenChoice_day2(i) == 1
                valChoWrong(i) = Mday1(idxCho_day2_R(i),3);
                valUncWrong(i) = Mday1(idxCho_day2_L(i),3);
            else
                valChoWrong(i) = NaN;
                valUncWrong(i) = NaN;
            end
        end
        difValWrong = abs(valChoWrong - valUncWrong); % value difference
        chMuncWrong = valChoWrong - valUncWrong; % value chosen - value unchosen
        
        %% movement onset and side
        movCho = Mday2(Mday2(:,3)==1,8);
        sidCho = Mday2(Mday2(:,3)==1,7);
        
        %% find indices of objects within the 120 used ones
        objs        = readtable([dir.behDat,fs,'Objects.csv']);
        foo         = logical(table2array(objs(:,7)));
        objsName    = table2cell(objs(foo,2)); clear foo objs
        
        for i = 1:length(idIma_day2)
            objIdx(i) = find(strcmp(objsNameAll{idIma_day2(i)},objsName));
        end
        
        %% remove nans if needed
        if excNan == 1
            objKeep = ~isnan(objVal) & ~isnan(objValWrong) & ~isnan(objCon) & ~isnan(objFam) & ~isnan(objPri);
            difKeep = ~isnan(difVal) & ~isnan(difValWrong) & ~isnan(valCho) & ~isnan(movCho) & ~isnan(sidCho);
            onsIma = onsIma(objKeep);
            objVal = objVal(objKeep);
            objValWrong = objValWrong(objKeep);
            objCon = objCon(objKeep);
            objFam = objFam(objKeep);
            objPri = objPri(objKeep);
            objSes = objSes(objKeep);
            objIdx = objIdx(objKeep);
            onsCho = onsCho(difKeep);
            difVal = difVal(difKeep);
            valCho = valCho(difKeep);
            valChoWrong = valChoWrong(difKeep);
            valUnc = valUnc(difKeep);
            chMunc = chMunc(difKeep);
            chMuncWrong = chMuncWrong(difKeep);
            movCho = movCho(difKeep);
            sidCho = sidCho(difKeep);
        end
        
        %% median splits
        objVal_pert = objVal + 0.00001*[1:length(objVal)]';
        medVal = median(objVal_pert);
        idxVal_H = objVal_pert > medVal;
        idxVal_L = objVal_pert < medVal;
        
        objCon_pert = objCon + 0.00001*[1:length(objVal)]';
        medCon = median(objCon_pert);
        idxCon_H = objCon_pert > medCon;
        idxCon_L = objCon_pert < medCon;
        
        objFam_pert = objFam + 0.00001*[1:length(objVal)]';
        medFam = median(objFam_pert);
        idxFam_H = objFam_pert > medFam;
        idxFam_L = objFam_pert < medFam;
        
        objPri_pert = objPri + 0.00001*[1:length(objVal)]';
        medPri = median(objPri_pert);
        idxPri_H = objPri_pert > medPri;
        idxPri_L = objPri_pert < medPri;
        
        difVal_pert = difVal + 0.00001*[1:length(difVal)]';
        medDV = median(difVal_pert);
        idxDV_H = difVal_pert > medDV;
        idxDV_L = difVal_pert < medDV;
        
        valCho_pert = valCho + 0.00001*[1:length(difVal)]';
        medCho = median(valCho_pert);
        idxCho_H = valCho_pert > medCho;
        idxCho_L = valCho_pert < medCho;
        
        chMunc_pert = chMunc + 0.00001*[1:length(difVal)]';
        medCmU = median(chMunc_pert);
        idxCmU_H = chMunc_pert > medCmU;
        idxCmU_L = chMunc_pert < medCmU;
        
        %% session type
        bData(subs(s)).sessType{r} = day2Order{r};
        
        %% create output struct
        
        %%%%%%%%%%%%%%%%%%%%%
        % imagination stuff %
        %%%%%%%%%%%%%%%%%%%%%
        
        bData(subs(s)).imagination(r).names = objSes;   % object names
        bData(subs(s)).imagination(r).objIdx = objIdx;  % real index (i.e. in the folder)
        bData(subs(s)).imagination(r).onset = onsIma;   % presentation onset
        
        % subjective evaluations
        bData(subs(s)).imagination(r).val = objVal;     % object value
        bData(subs(s)).imagination(r).valWrong = objValWrong;     % object value - wrong
        bData(subs(s)).imagination(r).con = objCon;     % object confidence in value
        bData(subs(s)).imagination(r).fam = objFam;     % object familiarity
        bData(subs(s)).imagination(r).pri = objPri;     % object price
        
        % median splits
        bData(subs(s)).imagination(r).k_val.low.onset = onsIma(idxVal_L);
        bData(subs(s)).imagination(r).k_val.high.onset = onsIma(idxVal_H);
        bData(subs(s)).imagination(r).k_con.low.onset = onsIma(idxCon_L);
        bData(subs(s)).imagination(r).k_con.high.onset = onsIma(idxCon_H);
        bData(subs(s)).imagination(r).k_fam.low.onset = onsIma(idxFam_L);
        bData(subs(s)).imagination(r).k_fam.high.onset = onsIma(idxFam_H);
        bData(subs(s)).imagination(r).k_pri.low.onset = onsIma(idxPri_L);
        bData(subs(s)).imagination(r).k_pri.high.onset = onsIma(idxPri_H);
        
        %%%%%%%%%%%%%%%%
        % choice stuff %
        %%%%%%%%%%%%%%%%
        
        bData(subs(s)).choice(r).onset = onsCho;
        bData(subs(s)).choice(r).valDiff = difVal;
        bData(subs(s)).choice(r).valCho = valCho;
        bData(subs(s)).choice(r).valChoWrong = valChoWrong;
        bData(subs(s)).choice(r).valUnc = valUnc;
        bData(subs(s)).choice(r).chMunc = chMunc;
        bData(subs(s)).choice(r).chMuncWrong = chMuncWrong;
        bData(subs(s)).choice(r).movOnset = movCho;
        bData(subs(s)).choice(r).movSide = sidCho;
        bData(subs(s)).choice(r).objPresented = choObjPresented;
        bData(subs(s)).choice(r).objChosen = chosenChoice_day2;
        
        % median splits
        bData(subs(s)).choice(r).k_valDiff.low.onset = onsCho(idxDV_L);
        bData(subs(s)).choice(r).k_valDiff.high.onset = onsCho(idxDV_H);
        bData(subs(s)).choice(r).k_valCho.low.onset = onsCho(idxCho_L);
        bData(subs(s)).choice(r).k_valCho.high.onset = onsCho(idxCho_H);
        bData(subs(s)).choice(r).k_chMunc.low.onset = onsCho(idxCmU_L);
        bData(subs(s)).choice(r).k_chMunc.high.onset = onsCho(idxCmU_H);
        
    end
end