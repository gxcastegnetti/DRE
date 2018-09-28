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
objs        = readtable([dir.beh,fs,'Objects.csv']);
objsNameAll = table2cell(objs(:,2)); clear foo objs

%% loop over subjects and sessions
for s = 1:length(subs)
    for r = 1:n_sess
        
        %% extract subject-specific task order
        if taskOrd(s) == 1
            day2Order = {'F','B','F','B'};
            day1Order = {'1','2','1','2'};
        elseif taskOrd(s) == 2
            day2Order = {'B','F','B','F'};
            day1Order = {'2','1','2','1'};
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
        difVal = abs(Mday1(idxCho_day2_L,3) - Mday1(idxCho_day2_R,3)); % value difference
        chMunc = valCho - valUnc; % value chosen - value unchosen
        
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
        objPri = Mday1F(idxFam_day2,4); % monetary value
        
        %% find indices of objects within the 120 used ones
        objs        = readtable([dir.beh,fs,'Objects.csv']);
        foo         = logical(table2array(objs(:,7)));
        objsName    = table2cell(objs(foo,2)); clear foo objs
        
        for i = 1:length(idIma_day2)
            objIdx(i) = find(strcmp(objsNameAll{idIma_day2(i)},objsName));
        end
        
        %% create output struct
        
        % remove NaNs
        if excNan == 1
            objKeep = ~isnan(objVal) & ~isnan(objCon) & ~isnan(objFam) & ~isnan(objPri);
            difKeep = ~isnan(difVal) & ~isnan(valCho) & ~isnan(movCho) & ~isnan(sidCho);
            onsIma = onsIma(objKeep);
            objVal = objVal(objKeep);
            objCon = objCon(objKeep);
            objFam = objFam(objKeep);
            objPri = objPri(objKeep);
            objSes = objSes(objKeep);
            objIdx = objIdx(objKeep);
            onsCho = onsCho(difKeep);
            difVal = difVal(difKeep);
            valCho = valCho(difKeep);
            valUnc = valUnc(difKeep);
            chMunc = chMunc(difKeep);
            movCho = movCho(difKeep);
            sidCho = sidCho(difKeep);
        end
        
        % median splits
        objVal = objVal + 0.00001*[1:length(objVal)]';
        medVal = median(objVal);
        idxVal_H = objVal > medVal;
        idxVal_L = objVal < medVal;
        
        objCon = objCon + 0.00001*[1:length(objVal)]';
        medCon = median(objCon);
        idxCon_H = objCon > medCon;
        idxCon_L = objCon < medCon;
        
        objFam = objFam + 0.00001*[1:length(objVal)]';
        medFam = median(objFam);
        idxFam_H = objFam > medFam;
        idxFam_L = objFam < medFam;
        
        objPri = objPri + 0.00001*[1:length(objVal)]';
        medPri = median(objPri);
        idxPri_H = objPri > medPri;
        idxPri_L = objPri < medPri;
        
        difVal = difVal + 0.00001*[1:length(difVal)]';
        medDV = median(difVal);
        idxDV_H = difVal > medDV;
        idxDV_L = difVal < medDV;
        
        valCho = valCho + 0.00001*[1:length(difVal)]';
        medCho = median(valCho);
        idxCho_H = valCho > medCho;
        idxCho_L = valCho < medCho;
        
        chMunc = chMunc + 0.00001*[1:length(difVal)]';
        medCmU = median(chMunc);
        idxCmU_H = chMunc > medCmU;
        idxCmU_L = chMunc < medCmU;
        
        % session type
        bData(subs(s)).sessType{r} = day2Order{r};
        
        %%%%%%%%%%%%%%%%%%%%%
        % imagination stuff %
        %%%%%%%%%%%%%%%%%%%%%
        
        bData(subs(s)).imagination(r).names = objSes;   % object names
        bData(subs(s)).imagination(r).objIdx = objIdx;  % real index (i.e. in the folder)
        bData(subs(s)).imagination(r).onset = onsIma;   % presentation onset
        
        % subjective evaluations
        bData(subs(s)).imagination(r).val = objVal;     % object value
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
        bData(subs(s)).choice(r).valUnc = valUnc;
        bData(subs(s)).choice(r).chMunc = chMunc;
        bData(subs(s)).choice(r).movOnset = movCho;
        bData(subs(s)).choice(r).movSide = sidCho;
        
        % median splits
        bData(subs(s)).choice(r).k_valDiff.low.onset = onsCho(idxDV_L);
        bData(subs(s)).choice(r).k_valDiff.high.onset = onsCho(idxDV_H);
        bData(subs(s)).choice(r).k_valCho.low.onset = onsCho(idxCho_L);
        bData(subs(s)).choice(r).k_valCho.high.onset = onsCho(idxCho_H);
        bData(subs(s)).choice(r).k_chMunc.low.onset = onsCho(idxCmU_L);
        bData(subs(s)).choice(r).k_chMunc.high.onset = onsCho(idxCmU_H);
        
    end
end