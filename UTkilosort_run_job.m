function status=UTkilosort_run_job(job_file)
% status=UTkilosort_run_job(job_file)
    %runs kilosort using input parameters located in job_file.
    %puts results in job.results_path and [job.results_path,'_after_automerge'] 
 global MULTICHANNEL_SORTING_PATH   
 if(isempty(MULTICHANNEL_SORTING_PATH))
     error('MULTICHANNEL_SORTING_PATH is not defined. This should be defined in BaphyConfigPath.')
 end
 addpath(genpath([MULTICHANNEL_SORTING_PATH, 'KiloSort'])) % path to kilosort folder
 
 do_write=true;
status=0;
job=load(job_file);
if exist([job.results_path filesep],'dir')
    [s,m,mid]=rmdir([job.results_path filesep]);
    if ~s, error(m), end
end
UTmkdir([job.results_path filesep])
[w,s]=unix(['chmod 777 ',job.results_path]);if w, error(s), end
slashes=strfind(job.results_path,filesep);
for sli=0:1
    [w,s]=unix(['chmod 777 ',job.results_path(1:slashes(end-sli))]);if w, warning(s), end
end
UTmkdir([job.results_path_temp filesep])
[w,s]=unix(['chmod 777 ',job.results_path_temp]);if w, error(s), end
for i=1:length(job.runs)
    LoadMFile([job.runs_root filesep job.runs{i}]);
    globalparams_{i}=globalparams;
    if strcmp(globalparams.HWparams.DAQSystem,'Open-Ephys')
        job.runinfo(i)=rawgetinfo([job.runs_root filesep job.runs{i}],globalparams);
        run_pref=strsep(job.runs{i},'_');
        run_pref=run_pref{1};
        checktgzevp=[job.runs_root filesep 'raw' filesep run_pref '.tgz'];
        tic;
        evpfile=evpmakelocal(checktgzevp,job.runinfo(i).evpv);
        slashes=strfind(evpfile,filesep);
        for sli=0:3
            [w,s]=unix(['chmod 777 ',evpfile(1:slashes(end-sli))]); if w, warning(s), end
        end
        fprintf('evpmakelocal for %s took %3.0fs. \n', job.runs{i},toc);
    else
        job.exptparams{i}=exptparams;
        job.globalparams(i)=globalparams;
        job.runinfo(i)=rawgetinfo(globalparams.evpfilename,globalparams);
        run_pref=strsep(job.runs{i},'_');
        run_pref=run_pref{1};
        checktgzevp=[job.runs_root filesep 'raw' filesep run_pref '.tgz'];
        evpfile=evpmakelocal(checktgzevp,job.runinfo(i).evpv);
    end
    if(isempty(evpfile))
       error('')
    end
    job.root{i}=evpfile;
    job.StartTime(i,:)=exptparams.StartTime;
    job.StartTime_re_Run1(i)=etime(job.StartTime(i,:),job.StartTime(1,:));
end
time_re_sound_target=.3;
min_offset_gap=.2;
switch job.runinfo(1).evpv
    case 5
        recording_start=[0 cumsum(job.nSamplesBlocks(1:end-1))]+1;
        trial_onsets_= job.trial_onsets_;
        runs_per_trial=cellfun('length',trial_onsets_);
        for i=1:length(job.runs)
            trial_onsets_{i}=trial_onsets_{i}+sum(job.nSamplesBlocks(1:i-1));
        end
    case 6 %6 means Open-Ophys format
        for i=1:length(job.runs)
            not_filesep=setdiff({'/','\'},filesep);not_filesep=not_filesep{1};
            globalparams_{i}.rawfilename=strrep(globalparams_{i}.rawfilename,not_filesep,filesep);
            
            [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
            recording_start(i)=EVtimestamps{i}(1);
            trial_onsets_{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-recording_start(i)+1)'; % in samples re run start
            
            silence_between_stimuli=exptparams.TrialObject.ReferenceHandle.PreStimSilence+exptparams.TrialObject.ReferenceHandle.PostStimSilence;
            if (silence_between_stimuli-time_re_sound_target)>=min_offset_gap
                time_re_sound=time_re_sound_target;
            else
                time_re_sound=silence_between_stimuli-min_offset_gap;
                if time_re_sound<=0
                    warning('very little time between stimuli!')
                    time_re_sound=.05;
                end
            end
            drift_track_windows{i}=repmat(trial_onsets_{i}',1,2)-repmat(EVinfo(i).header.sampleRate*([time_re_sound,0]-exptparams.TrialObject.ReferenceHandle.PreStimSilence),length(trial_onsets_{i}),1);
            %trial_offsets=EVtimestamps(EVinfo.eventId==0&EVinfo.eventType==3)-recording_start+1; % in samples re recording start
            runs_per_trial(i)=length(trial_onsets_{i});
        end
end
trial_onsets=[trial_onsets_{:}];
job.drift_track_windows=drift_track_windows;
UTmkdir(job.fbinary)
if job.find_drift_correction
    job2=job;
    job2.driftCorrectionMode='None';
    if strcmp(job.driftCorrectionMode,'BeforeFiltering')
        job2.return_after_finding_drift_correction=true;
    end
    if isfield(job,'fbinary_uncorrected')
        job2.fbinary=job.fbinary_uncorrected;
    end
    [rez, DATA, uproj,job2] = preprocessData(job2,do_write,drift_track_windows);
    job.find_drift_correction=false;
    clear DATA rez
    [rez, DATA, uproj,job] = preprocessData(job,do_write,drift_track_windows);    
    job.find_drift_correction=true;
else
    if job.keep_fproc
        job.ForceMaxRAMforDat = 0;
    end
    [rez, DATA, uproj,job] = preprocessData(job,do_write,drift_track_windows);
end

%% save block (different concatenated files) sizes and start times
blocksizes=job.nSamplesBlocks;
blockstarts=job.StartTime_re_Run1*job.fs;
%actual saving below because rezToPhy deletes files
%% fit model
if ~job.GPU
    rez.ops.initialize=0;
end
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% save python results file for Phy
rezToPhy(rez, job.results_path);
save([job.results_path_temp filesep 'rez.mat'],'-Struct','rez', '-v7.3');
fprintf('Kilosort took %2.2f seconds\n', toc)

% now fire up Phy and check these results. 

%% AUTO MERGES 
rez = merge_posthoc2(rez);

%Fix merge problems:
%the automerge function merges clusters with less than a certain number
%of spikes and sets their IDs to 0. This ends up misaligning all the \
%similarity values. This code sets this cluster's ids to an empty slot,
%and shifts the other ids down by 1 to start at 0.
if(~isequal(rez.st3(:,2),rez.st3(:,5)))
    %ids are shifted off by one in the merged version
    
    %rez.st3(:,5) contains the cluster identities for each spike
    %shift indexes down 1
    rez.st3(:,5)=rez.st3(:,5)-1;
    
    %find empty slots
    un=unique(rez.st3(:,5));
    empties=setdiff(1:max(un),un);
    if isempty(empties), error('why empty?'); end
    
    %set cluster that did have id=0 to an empty slot
    rez.st3(rez.st3(:,5)==-1,5)=empties(1);
end

% save python results file for Phy
UTmkdir([job.results_path,'_after_automerge',filesep]);
[w,s]=unix(['chmod 777 ',job.results_path,'_after_automerge']);
if w, error(s), end
rezToPhy(rez, [job.results_path,'_after_automerge']);

%% save and clean up
% save matlab results file for future use (although you should really only be using the manually validated spike_clusters.npy file)
rez_after.st3=rez.st3;
save([job.results_path_temp filesep 'rez_after_automerge.mat'],'-Struct','rez_after', '-v7.3');

%save blocksizes and blockstarts
writeNPY(blocksizes,[job.results_path,'_after_automerge',filesep,'blocksizes.npy']);
writeNPY(blockstarts,[job.results_path,'_after_automerge',filesep,'blockstarts.npy']);
writeNPY(blocksizes,[job.results_path,filesep,'blocksizes.npy']);
writeNPY(blockstarts,[job.results_path,filesep,'blockstarts.npy']);

%give everyone permission to do everything. Anarchy, whooo!
[w,s]=unix(['chmod -R 777 ',job.results_path]);if w, error(s), end
[w,s]=unix(['chmod -R 777 ',job.results_path,'_after_automerge']);if w, error(s), end
[w,s]=unix(['chmod -R 777 ',job.results_path_temp]);if w, error(s), end

% remove temporary file
delete(job.fproc);
job.status=1;
job.kcoords=rez.ops.kcoords;
job.nt0min=rez.ops.nt0min;
job.min_spike_perc_keep=rez.ops.min_spike_perc_keep;
save(job_file,'-Struct','job')
status=1;
end