function status=UTkilosort2_run_job(job_file)
% status=UTkilosort_run_job(job_file)
    %runs kilosort using input parameters located in job_file.
    %puts results in job.results_path and [job.results_path,'_after_automerge'] 
global MULTICHANNEL_SORTING_PATH BAPHYDATAROOT LOCAL_DATA_ROOT
if(isempty(MULTICHANNEL_SORTING_PATH))
 error('MULTICHANNEL_SORTING_PATH is not defined. This should be defined in BaphyConfigPath.')
end
do_write=true;

%% Check available memory
dev=gpuDevice(1);
freeGB = dev.AvailableMemory/1024/1024/1024;
totalGB = dev.TotalMemory/1024/1024/1024;
fprintf('%0.2g/%0.2g GB free on GPU.\n', freeGB,totalGB);
if freeGB < 3
    [~, name] = system('hostname'); name(end)=[];
    status = mysql(['UPDATE tComputer SET maxGPU_jobs=0 WHERE name=''',name,''';']);
    error(['Not enough memory free, killing job, setting maxGPU_jobs to 0 for ' name ' to stop more sorting jobs from going here.',...
        ' Run the following in mysql to turn it back on: UPDATE tComputer SET maxGPU_jobs=1 WHERE name=''',name,''';'])
end

%% Create local copy of data
status=0;
job=load(job_file);
if exist([job.results_path filesep],'dir')
    [s,m,mid]=rmdir([job.results_path filesep], 's');
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



UTmkdir(job.results_path);
[w,s]=unix(['chmod 777 ', job.results_path]); if w, error(s), end


% creates a local binary file and "swap" to expedite procecing. 
% Stores e the server binary location for later recovery.
job.fbinary_server = job.fbinary;
job.fbinary = strrep(job.fbinary, BAPHYDATAROOT, LOCAL_DATA_ROOT);
job.fproc=[fileparts(job.fbinary) filesep 'temp_wh.dat'];% "swap" or RAM replacement

UTmkdir(fileparts(job.fbinary)); % this should create folder for fproc too
UTmkdir(fileparts(job.fbinary_server));

% anarchy permissions.
[w,s]=unix(['chmod 777 ',fileparts(job.fbinary)]); if w, warning(s), end
[w,s]=unix(['chmod 777 ',fileparts(job.fbinary_server)]); if w, warning(s), end

% related with Luke drift correction, does it intervene in something else?
% todo: delete
%if job.keep_fproc
%    job.ForceMaxRAMforDat = 0;
%end

%% Kilosort2 master start. MLE.

addpath(genpath([MULTICHANNEL_SORTING_PATH, 'KiloSort2'])) % path to kilosort folder
addpath(genpath([MULTICHANNEL_SORTING_PATH, 'open-ephys-matlab-tools'])) % perhaps this has all what we need! MLE
rmpath(genpath([MULTICHANNEL_SORTING_PATH, 'KiloSort2' filesep 'CUDA'])) % path to kilosort folder
if ~verLessThan('matlab', '9.8')%2020a
    addpath([MULTICHANNEL_SORTING_PATH, 'KiloSort2' filesep 'CUDA' filesep 'MATLAB2020'])
%elseif ~verLessThan('matlab', '9.3')%2017a
else
    addpath([MULTICHANNEL_SORTING_PATH, 'KiloSort2' filesep 'CUDA' filesep 'MATLAB2017'])
    %Might not work for versions < 2017a, need to recompile (re-run mexGPUall)
end
    
addpath([MULTICHANNEL_SORTING_PATH, 'npy-matlab']) % path to npy-matlab

%% LBHB preprocessing

switch job.datatype
    case 'Open-Ephys'
        % old OPE hav speciale format, new OEP are binary files
        if strcmp(job.runinfo(1).datatype, 'binary')
            job = concatenateOpenEphysBinary(job);
        else
            job = convertOpenEphysToRawBInary(job);
        end
        % transforms the channel map form file to struct and sets into 1
        % indexing, savese the channelMap directory for cleaning up
        job.chanMapDir = job.chanMap;
        job.chanMap = load(job.chanMap);
        job.chanMap.chanMap = job.chanMap_KiloRaw;
        
    case 'MANTA'
        job = convertMANTAToRawBinary(job,do_write);  % convert data, only for MANTA
        job.chanMap = chanMap;   
end

%% Get trial onset times, 
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
            
            if strcmp(job.runinfo(i).datatype, 'binary')
                %Get trial onsets from events in continuous folder
                EVdata = load_open_ephys_binary(job.runinfo(i).json_file,'events',job.runinfo(i).event_ind_TTL);
                sampleRate(i) = EVdata.Header.sample_rate;
                job.recording_start(i) = EVdata.Timestamps(1);
                
                % onsets in samples from the recording start
                trial_onsets_{i}=(EVdata.Timestamps(EVdata.Data==1) - ...
                                  job.recording_start(i)+1)' + ...
                                  sum(job.nSamplesBlocks(1:i-1)); 
               
            else
                [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
                sampleRate(i)=EVinfo(i).header.sampleRate;
                %recording_startO(i)=EVtimestamps{i}(1);

                data_files=dir([job.root{i},filesep,'*.continuous']);
                [~, EVinfo_(i)] = load_open_ephys_data_faster([job.root{i},filesep,data_files(1).name],0);
                job.recording_start(i)=EVinfo_(i).ts(1);

                trial_onsets_{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-job.recording_start(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start
                %trial_onsets_O{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-recording_startO(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re run start

            end
            
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
            runs_per_trial(i)=length(trial_onsets_{i});
        end
end
trial_onsets=[trial_onsets_{:}];
job.trial_onsets_ = trial_onsets_;
job.trial_onsets = trial_onsets;
job.runs_per_trial = runs_per_trial;
job.sampleRate = sampleRate;


%% Runs the photoelectric artifact corrections

% check if need to remove laser artifacts
remove_laser_artifact_sec = getparm(job, 'remove_laser_artifact_sec', 0);
interp_laser_sec = getparm(job, 'interp_laser_sec', 0);
if remove_laser_artifact_sec>0
   
    fprintf('correcting photoelectric artifact\n')
    fprintf('removing %.3f sec around on/off\n', ...
            remove_laser_artifact_sec);
    if interp_laser_sec>0
      fprintf('interpolating %.5f sec around on/off\n', interp_laser_sec);
    end
    
    tic
    % need to get experiments with laser stim and  events to figure out 
    % when laser was turned on/off, so supply baphy events file
    ncorr = 0;
    [tmp pid] = system('pgrep MATLAB');
    for j=1:length(job.runs)
        
        
        [tmp mem_usage] = system(['cat /proc/' strtrim(pid) '/status | grep VmData']);
        fprintf("prior loading %i MB\n", round(str2num(strtrim(extractAfter(extractBefore(mem_usage, ' kB'), ':'))) / 1000));

        
        parmfile = fullfile(job.runs_root, job.runs{j});
        parmout = LoadMFile(parmfile);
        exptparams = parmout.exptparams;
        exptevents = parmout.exptevents;
        
        if ~strcmp(exptparams.TrialObjectClass,'RefTarOpt')
            % skips no opto experiments
            continue
        else
            fprintf('block %d\n', j);
            ncorr = ncorr + 1;
        end 
               
        spikefs = job.fs;
        
        % reads from the big binary file, the chunk corresponding to the
        % experiment block
        fid = fopen(job.fbinary, 'r+');  
        offset = job.Nchan * sum(job.nSamplesBlocks(1:j-1)) * 2;
        fseek(fid, offset, 'bof');
        all_samples = fread(fid, [job.Nchan, job.nSamplesBlocks(j)], 'int16');
        

        [tmp mem_usage] = system(['cat /proc/' strtrim(pid) '/status | grep VmData']);
        fprintf("after loading %i MB\n", round(str2num(strtrim(extractAfter(extractBefore(mem_usage, ' kB'), ':'))) / 1000));

        
        for chan = 1:size(all_samples,1)
          tc = all_samples(chan,:)';
          tc_out = remove_opto_artifacts(tc,...
             job.trial_onsets_{j} -sum(job.nSamplesBlocks(1:j-1)),...
             spikefs, exptparams, exptevents, ...
             remove_laser_artifact_sec, interp_laser_sec);
          all_samples(chan,:) = int16(tc_out)';
        end
        
        % writes back to the same file location
        fseek(fid, offset, 'bof');
        fwrite(fid, all_samples, 'int16');
        fclose(fid);
        clear('all_samples');

        
        [tmp mem_usage] = system(['cat /proc/' strtrim(pid) '/status | grep VmData']);
        fprintf("after clearing %i MB\n", round(str2num(strtrim(extractAfter(extractBefore(mem_usage, ' kB'), ':'))) / 1000));

    end
    if ncorr > 0
        fprintf('done correcting photoelectric artifact, ')
    else
        fprintf('no RefTarOpt experiments found')
    end
    toc
    fprintf('\n')
end

%%% End LBHB special code

%% LBHB save block (different concatenated files) sizes and start times
blocksizes=job.nSamplesBlocks;
blockstarts=job.StartTime_re_Run1*job.chanMap.fs;
% actual saving below because rezToPhy deletes files

% the binary file is in this folder
rootZ = job.results_path_temp;

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
% fs = dir(fullfile(rootZ, 'chan*.mat'));
% if ~isempty(fs) && ~isfield(job, 'chanMap')
%     job.chanMap = fullfile(rootZ, fs(1).name);
% end

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(job);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
UTkilosort2_plot_drift(rez,1); % plot drift and save in results file
%save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');


% main tracking and template matching algorithm
try
    rez = learnAndSolve8b(rez);
catch err
    dev=gpuDevice(1);
    freeGB = dev.AvailableMemory/1024/1024/1024;
    totalGB = dev.TotalMemory/1024/1024/1024;
    fprintf('%0.2g/%0.2g GB free on GPU.\n', freeGB,totalGB);
    fprintf(err)
    rethrow(err)
end
% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))


%% LBHB save results to phy (MLE)
% copy binary to server and restablis pointer to files
[w,s]=unix(['cp ', job.fbinary, ' ', job.fbinary_server]); if w, warning(s), end
job.fbinary = job.fbinary_server;
job = rmfield(job, 'fbinary_server');

fprintf('Saving results to Phy in %s  \n', job.results_path)
rezToPhy(rez, job.results_path);
% save([job.results_path_temp filesep 'rez.mat'],'-Struct','rez', '-v7.3');
% fprintf('Kilosort took %2.2f seconds\n', toc)

%% if you want to save the results to a Matlab file... 
if 0
    % discard features in final rez file (too slow to save)
    rez.cProj = [];
    rez.cProjPC = [];

    % transforms GPU arrays into normal arrays
    fn=fieldnames(rez);
    for i=1:length(fn)
        if isa(rez.(fn{i}),'gpuArray')
            rez.(fn{i})=gather_try(rez.(fn{i}));
        end
    end
    fn=fieldnames(rez.ops);
    for i=1:length(fn)
        if isa(rez.ops.(fn{i}),'gpuArray')
            rez.ops.(fn{i})=gather_try(rez.ops.(fn{i}));
        end
    end

    % save final results as rez2
    fname = [job.results_path_temp filesep 'rez.mat'];
    fprintf('Saving final results in %s  \n', fname);
    save(fname,'-Struct', 'rez', '-v7.3');
end
% Kilosort2 master end. MLE
%% save and clean up

% save blocksizes and blockstarts
writeNPY(blocksizes,[job.results_path,filesep,'blocksizes.npy']);
writeNPY(blockstarts,[job.results_path,filesep,'blockstarts.npy']);

%give everyone permission to do everything. Anarchy, whooo!
[w,s]=unix(['chmod -R 777 ',job.results_path]);if w, error(s), end
% no automerge needed with KS2
%[w,s]=unix(['chmod -R 777 ',job.results_path,'_after_automerge']);if w, error(s), end
[w,s]=unix(['chmod -R 777 ',job.results_path_temp]);if w, error(s), end

% re-establish chanMap to default i.e. path.
job.chanMap = job.chanMapDir;
job = rmfield(job, 'chanMapDir');

% remove temporary file
delete(job.fproc);

job.status=1;
job.kcoords=rez.ops.kcoords;
job.nt0min=rez.ops.nt0min;
%job.min_spike_perc_keep=rez.ops.min_spike_perc_keep;
save(job_file,'-Struct','job')
end