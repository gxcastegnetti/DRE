function physioRegressors(dir,subs)
%% creates physiological regressors
% ~~~
% GX Castegnetti --- 2018

fs = filesep;
n_sess = 4;

nslices   = 48;
ndummies  = 5;  % Number of scans excluded from analysis
TR_sl     = 0.070;       % Slice TR in secs
TRms      = TR_sl*1e3;     % As above
nsessions = 1; % Number of scanning sessions in the file
slicenum  = 24; % Slice number to time-lock regressors to
sliceord  = 'ascending'; % Ascending by default.
slicenum  = get_slicenum(slicenum,nslices,sliceord);

for s = 1:length(subs)
    
    disp(['Creating physiological regressors for sub#', num2str(subs(s),'%03d'),'...']);
    
    % director with physiological data
    dirPhy = [dir.dre,'data',fs,'physio',fs,num2str(subs(s),'%03d')];
    
    % get Spike2 outputs
    for r = 1:n_sess
        
        % possible names of the file
        smrFileF = [dirPhy,fs,'SF',num2str(subs(s),'%03d'),'_S',num2str(r),'F.smr'];
        smrFileB = [dirPhy,fs,'SF',num2str(subs(s),'%03d'),'_S',num2str(r),'B.smr'];
        
        % take the correct one
        if exist(smrFileF,'file')
            smrFile = smrFileF;
        elseif exist(smrFileB,'file')
            smrFile = smrFileB;
        else
            warning(['sub#',num2str(subs(s),'%03d'),', S',num2str(r),' ---> None of the files (F,B) exists.'])
            continue
        end
        
        % set (and check) channels
        show_channels(smrFile);
        scanner_channel    = 1;
        cardiacTTL_channel = 2;
        cardiacQRS_channel = [];
        resp_channel       = 4;
        check_channels(smrFile,scanner_channel,cardiacTTL_channel,cardiacQRS_channel,resp_channel);
        
        % Call the main routine for calculating physio regressors
        % NB - currently the cardiacqrs calculation is disabled.
        try
            [cardiac,~,respire,rvt] = make_physio_regressors(smrFile,nslices,ndummies,TR_sl,...
                slicenum,nsessions,scanner_channel,cardiacTTL_channel,cardiacQRS_channel,resp_channel);
        catch
            continue
        end
        
        for sessnum = 1:nsessions
            R = [];
            if ~isempty(cardiac{sessnum})
                cardiac_sess = cardiac{sessnum};
                filename = sprintf('%s_cardiac_S%d',spm_str_manip(smrFile,'r'),sessnum);
                save(filename, 'cardiac_sess');
                R=cat(2,R,cardiac{sessnum}(:,1:6));
            end
            if ~isempty(respire) && ~isempty(respire{sessnum})
                respire_sess = respire{sessnum};
                filename = sprintf('%s_respire_S%d',spm_str_manip(smrFile,'r'),sessnum);
                save(filename, 'respire_sess');
                R=cat(2,R,respire{sessnum}(:,1:6));
            end
            if ~isempty(rvt) && ~isempty(rvt{sessnum})
                rvt_sess = rvt{sessnum};
                filename = sprintf('%s_rvt_S%d',spm_str_manip(smrFile,'r'),sessnum);
                save(filename,'rvt_sess');
                R=cat(2,R,rvt{sessnum}(:,1:size(rvt{sessnum},2)));
            end
            nfiles=size(R,1);
            
            % Save R for all physio only
            if nfiles>0
                oR=R;
                Rname = sprintf('%s_R_S%d',spm_str_manip(smrFile,'r'),sessnum);
                R=R-repmat(mean(R),nfiles,1);
                save(Rname, 'R');
            end
            % if required, concatenate with motion parameters
            
        end
    end
end
