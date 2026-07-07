function status=UTkilosort3_run_job(job_file)
% status=UTkilosort_run_job(job_file)
    %runs kilosort using input parameters located in job_file.
    %puts results in job.results_path and [job.results_path,'_after_automerge'] 
global MULTICHANNEL_SORTING_PATH
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

status=0;
job=load(job_file);


%% Create local copy of data

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

%% creates all the local relevant directories, saving the location of server counterpart

% todo figure out a more streamline way of passing these paths

if isfield(job, 'keep_local') && job.keep_local
    
    % if reruning a job, server locations probably already exist and we
    % dont want o override them with the local directions
    
    if ~isfield(job, 'fbinary_server') && ~isfield(job, 'results_path_server')...
            && ~isfield(job, 'results_path_temp_server')
        
        job.fbinary_server = job.fbinary;
        job.results_path_server = job.results_path;
        job.results_path_temp_server = job.results_path_temp;


        job.fbinary =  strrep(job.fbinary, BAPHYDATAROOT, '/KiloSort/');
        job.results_path = strrep(job.results_path, BAPHYDATAROOT, '/KiloSort/');
        job.results_path_temp = strrep(job.results_path_temp, BAPHYDATAROOT, '/KiloSort/');
    end
else
    job.keep_local=0;
end

job.fproc = [job.results_path_temp filesep 'temp_wh.dat'];

fprintf("creating local files at %s \n", fileparts(job.fbinary));
UTmkdir(fileparts(job.fbinary)); % this should create folder for fproc too

if exist([job.results_path filesep],'dir')
    [s,m,mid]=rmdir([job.results_path filesep], 's');
    if ~s, error(m), end
end
UTmkdir([job.results_path filesep])

% anarchy permissions.
[w,s]=unix(['chmod 777 ',fileparts(job.fbinary)]); if w, warning(s), end
[w,s]=unix(['chmod 777 ',job.results_path]);if w, error(s), end

if job.keep_local
    UTmkdir(fileparts(job.fbinary_server));
    [w,s]=unix(['chmod 777 ',fileparts(job.fbinary_server)]); if w, warning(s), end
end


%% Kilosort3 master paths. MLE.

addpath(genpath([MULTICHANNEL_SORTING_PATH, 'KiloSort3'])) % path to kilosort folder
addpath([MULTICHANNEL_SORTING_PATH, 'npy-matlab']) % path to npy-matlab

%% LBHB preprocessing

switch job.datatype
    case 'Open-Ephys'
        % old OPE hav special format, new OEP are binary files
        if strcmp(job.runinfo(1).datatype, 'binary')
            job = concatenateOpenEphysBinary(job);
        elseif strcmp(job.runinfo(1).datatype, 'Binary')
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
                sampleRate(i)=EVdata.Header.sample_rate;
                
                loc_json = strrep(job.runinfo(i).json_file, BAPHYDATAROOT, LOCAL_DATA_ROOT);
                continuous_data = load_open_ephys_binary(loc_json,'continuous',job.runinfo(i).data_ind,'mmap');
                job.recording_start(i) = continuous_data.Timestamps(1);
                
                % onsets in samples from the recording start
                trial_onsets_{i}=(EVdata.Timestamps(EVdata.Data==1) - ...
                                  job.recording_start(i)+1)' + ...
                                  sum(job.nSamplesBlocks(1:i-1)); 
               
            elseif strcmp(job.runinfo(i).datatype, 'Binary')

                spike_rec = job.runinfo(i).spike_rec;
                sampleRate(i)=job.runinfo(i).spikefs;
                trial_onsets = job.runinfo(i).recordings(spike_rec).trial_onsets;
                trial_offsets = job.runinfo(i).recordings(spike_rec).trial_offsets;
                
                cont = job.runinfo(i).session.recordNodes{spike_rec}.recordings{1}.continuous;
                ev = job.runinfo(i).session.recordNodes{spike_rec}.recordings{1}.ttlEvents;
                event_keys = ev.keys();
                k = event_keys{find(contains(event_keys,'AP') & contains(event_keys,'Neuropix'),1)};
                cont(k).metadata.startTimestamp
                
                % onsets in samples from the recording start
                trial_onsets_{i}=(trial_onsets - cont(k).metadata.startTimestamp+1)' + ...
                                  sum(job.nSamplesBlocks(1:i-1)); 
               
            else
                x = dir([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
                if isempty(x)
                    x=dir([globalparams_{i}.rawfilename,filesep,'*',filesep,'all_channels.events']);
                end
                OEfolder=x(1).folder;
                [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([OEfolder,filesep,'all_channels.events']);
                
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

%% LBHB save block (different concatenated files) sizes and start times
blocksizes=job.nSamplesBlocks;
blockstarts=job.StartTime_re_Run1*job.chanMap.fs;
% actual saving below because rezToPhy deletes files

% the binary file is in this folder
rootZ = job.results_path_temp;

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)
disp(job)
% % main parameter changes from Kilosort2 to v2.5
% job.sig        = 20;  % spatial smoothness constant for registration
% job.fshigh     = 300; % high-pass more aggresively
% job.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
% 
% % main parameter changes from Kilosort2.5 to v3.0
% ops.Th       = [9 9];

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if 0 %~isempty(fs) && ~isfield(job, 'chanMap')
    job.chanMap = fullfile(rootZ, fs(1).name);
end

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(job);

% NEW STEP TO DO DATA REGISTRATION
rez = datashift2(rez, 1); % last input is for shifting data
set(gcf,'Color','w');
saveas(gcf,[job.results_path filesep 'drift_map'],'fig')
set(findall(gcf,'LineStyle','--'),'Visible','off')
export_fig([job.results_path filesep 'drift_map.png'],'-png','-a1','-m3','-q100')
if getOr(job, 'recompute_drift_map_after_shifting', 1)
    rez_recalc=rmfield(rez,{'wTEMP','wPCA','iC','dist','iorig','dshift','st0','F','F0','F0m'});
    rez_recalc.ops=rmfield(rez_recalc.ops,{'dmin','yup','dminx','xup'});
    rez_recalc = datashift2(rez_recalc,0);
    fns={'dshift','F','F0'};
    for i=1:length(fns)
        rez.post_shift.(fns{i})=rez_recalc.(fns{i});
    end
    title('Post-shift drift map')
    set(gcf,'Color','w');
    saveas(gcf,[job.results_path filesep 'drift_map_post_shift'],'fig')
    set(findall(gcf,'LineStyle','--'),'Visible','off');
    export_fig([job.results_path filesep 'drift_map_post_shift.png'],'-png','-a1','-m3','-q100')
end

[rez, st3, tF]     = extract_spikes(rez);

rez                = template_learning(rez, tF, st3);

[rez, st3, tF]     = trackAndSort(rez);

rez                = final_clustering(rez, tF, st3);

% final merges
rez = find_merges(rez, 1);


%% LBHB save results to phy (MLE)
fprintf('Saving results to Phy in %s  \n', job.results_path)
rezToPhy2(rez, job.results_path);
% save([job.results_path_temp filesep 'rez.mat'],'-Struct','rez', '-v7.3');
% fprintf('Kilosort took %2.2f seconds\n', toc)

%% if you want to save the results to a Matlab file... 

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

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
%delete(job.fproc);

job.status=1;
job.kcoords=rez.ops.kcoords;
job.nt0min=rez.ops.nt0min;
%job.min_spike_perc_keep=rez.ops.min_spike_perc_keep;
save(job_file,'-Struct','job')
end