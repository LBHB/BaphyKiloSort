function job=UTkilosort2_fix_trial_onsets(job_file, do_save)

job=load(job_file);
if strcmp(job.datatype,'MANTA')
    fprintf('\nMANTA job, not doing anything')
    return
end
job = rmfield(job,'runinfo');
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
    fprintf(' %d',job.runs{i})
    bad_trials{i}=[];
    timestamp_skip_times{i}=[];
    timestamp_skip_trials{i}=[];
    timestamp_skip_gap_lengths{i}=[];
    
    not_filesep=setdiff({'/','\'},filesep);not_filesep=not_filesep{1};
    globalparams_{i}.rawfilename=strrep(globalparams_{i}.rawfilename,not_filesep,filesep);
    
    [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
    sampleRate(i)=EVinfo(i).header.sampleRate;
    recording_startO(i)=EVtimestamps{i}(1);
    
    data_files=dir([job.root{i},filesep,'*.continuous']);
    [~, EVinfo_(i),timestamps] = load_open_ephys_data_faster([job.root{i},filesep,data_files(1).name],0);
    job.recording_start(i)=EVinfo_(i).ts(1);
    
    these_trial_onsets = EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3);
    for j=1:length(these_trial_onsets)
        %Convert from time to index in data vector
        onset_ind = find(these_trial_onsets(j)==timestamps);
        bad_onset=length(onset_ind)~=1;
        if bad_onset
            [~,onset_ind]=min(abs(these_trial_onsets(j)-timestamps));
            closest_diff = timestamps(onset_ind)-these_trial_onsets(j);
            warning('Timestamp matching onset of trial %d (by TTL pulse) not found, closest timestamp is %d samples re onset.',j,closest_diff)
            bad_trials{i} = [bad_trials{i} j];
            timestamp_skip_trials{i} = [timestamp_skip_trials{i} j];
            timestamp_skip_times{i} = [timestamp_skip_times{i} 0];
            timestamp_skip_gap_lengths{i} = [timestamp_skip_gap_lengths{i} closest_diff];
        end
        these_trial_onsets(j) = onset_ind;
        if j< length(these_trial_onsets)
            offset_ind = find(these_trial_onsets(j+1)==timestamps)-1;
            dft=diff(timestamps(onset_ind:offset_ind));
        else
            dft=diff(timestamps(onset_ind:end));
        end
        unique_diffs=unique(dft);
        if ~isequal(unique_diffs,1)
            if ~bad_onset
                %Already marked above if bad onset
                bad_trials{i} = [bad_trials{i} j];
            end
            unique_diffs(unique_diffs==1)=[];
            for k=1:length(unique_diffs)
                skip_times = find(dft==unique_diffs(k));
                timestamp_skip_trials{i} = [timestamp_skip_trials{i} j];
                timestamp_skip_times{i} = [timestamp_skip_times{i} skip_times];
                timestamp_skip_gap_lengths{i} = [timestamp_skip_gap_lengths{i} repmat(unique_diffs(k),1,length(skip_times))];
            end
        end
    end
    trial_onsets_{i} = these_trial_onsets' + sum(job.nSamplesBlocks(1:i-1));
    %trial_onsets_{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-job.recording_start(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start
    %trial_onsets_O{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-recording_startO(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start
    
%     silence_between_stimuli=exptparams.TrialObject.ReferenceHandle.PreStimSilence+exptparams.TrialObject.ReferenceHandle.PostStimSilence;
%     if (silence_between_stimuli-time_re_sound_target)>=min_offset_gap
%         time_re_sound=time_re_sound_target;
%     else
%         time_re_sound=silence_between_stimuli-min_offset_gap;
%         if time_re_sound<=0
%             warning('very little time between stimuli!')
%             time_re_sound=.05;
%         end
%     end
    %trial_offsets=EVtimestamps(EVinfo.eventId==0&EVinfo.eventType==3)-recording_start+1; % in samples re recording start
    runs_per_trial(i)=length(trial_onsets_{i});
end

trial_onsets=[trial_onsets_{:}];
job.trial_onsets_ = trial_onsets_;
job.trial_onsets = trial_onsets;
job.runs_per_trial = runs_per_trial;
job.sampleRate = sampleRate;
job.trial_onsets_calculated_with_timestamps=1;
job.bad_trials=bad_trials;
job.has_timestamp_skips = any(~cellfun(@isempty,bad_trials));
job.timestamp_skip_times=timestamp_skip_times;
job.timestamp_skip_trials=timestamp_skip_trials;
job.timestamp_skip_gap_lengths=timestamp_skip_gap_lengths;

if do_save==2
    if job.has_timestamp_skips
        %Make a backup copy
        fprintf('\nCopying %s to %s and saving new job struct in %s',job_file,strrep(job_file,'.mat','_bck.mat'),job_file)
        [status,msg]=copyfile(job_file,strrep(job_file,'.mat','_bck.mat'));
        if status==0
            error(msg)
        end
        save(job_file,'-Struct','job')
    end
elseif do_save==1
    fprintf('\nSaving new job struct in %s',job_file)
    save(job_file,'-Struct','job')
end