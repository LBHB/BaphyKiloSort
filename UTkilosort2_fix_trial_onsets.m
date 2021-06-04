function job=UTkilosort2_fix_trial_onsets(job_file)

job=load(job_file);

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

for i=1:length(job.runs)
    not_filesep=setdiff({'/','\'},filesep);not_filesep=not_filesep{1};
    globalparams_{i}.rawfilename=strrep(globalparams_{i}.rawfilename,not_filesep,filesep);
    
    [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
    sampleRate(i)=EVinfo(i).header.sampleRate;
    recording_startO(i)=EVtimestamps{i}(1);
    
    data_files=dir([job.root{i},filesep,'*.continuous']);
    [~, EVinfo_(i)] = load_open_ephys_data_faster([job.root{i},filesep,data_files(1).name],0);
    job.recording_start(i)=EVinfo_(i).ts(1);
    
    trial_onsets_{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-job.recording_start(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start
    trial_onsets_O{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-recording_startO(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start
    
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
    %trial_offsets=EVtimestamps(EVinfo.eventId==0&EVinfo.eventType==3)-recording_start+1; % in samples re recording start
    runs_per_trial(i)=length(trial_onsets_{i});
end

trial_onsets=[trial_onsets_{:}];
job.trial_onsets_ = trial_onsets_;
job.trial_onsets = trial_onsets;
job.runs_per_trial = runs_per_trial;
job.sampleRate = sampleRate;

save(job_file,'-Struct','job')

end