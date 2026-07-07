function kilosort_run_queue_job(parmfile)

%script for creating KiloSort jobs and adding them to the jobs folder
global BAPHYDATAROOT SERVER_PATH
if(isempty(BAPHYDATAROOT))
    baphy_set_path
end

[pp,bb,ee]=fileparts(parmfile);
if strcmpi(ee,'.json')
   fid = fopen(parmfile); 
   raw = fread(fid,inf); 
   str = char(raw'); 
   fclose(fid); 
   options = jsondecode(str);
else
   options = load(parmfile); 
end

options.animal=getparm(options,'animal','Tartufo');
options.site=getparm(options,'site','TAR017b');
options.probe=getparm(options,'probe','');
options.run_nums=getparm(options,'run_nums',[10]); %can be 'all'
options.Nfilt=getparm(options,'Nfilt',0);
options.common_rejection_mode=getparm(options,'common_rejection_mode','median');
options.KSversion=getparm(options,'KSversion',2);
options.skip_channels=getparm(options,'skip_channels','');
options.rm_noise=getparm(options,'rm_noise',0);
% if specified, options.channel_map will force use of a specific channel
% map if not specified correctly during data acquisition
% options.channel_map = 

% subtract average signal around each laser onset and offset? if non-zero,
% subtract over this many seconds following laser onset (suggest 0.1 sec?)
% this work is handled by baphy.Utiltities.Kilosort.remove_opto_artifacts
% called from
% multichannelsorting.Kilosort2.preprocessing.convertOpenEphysToRawBinary
options.remove_laser_artifact_sec=getparm(options,'remove_laser_artifact_sec',0);
% interpolate this many bins after subtraction (starting 5 bins before 
% laser onset--to avoid jitter issues) to really make sure no artifact 
% survived (typically for brief/bad noise, that lasts just ~0.0015 sec)
options.interp_laser_sec=getparm(options,'interp_laser_sec',0.0015);

KSpaths=UTkilosort_paths();


% All this code is to fill in job.runs_root and job.runs
if options.KSversion==2
    job=UTkilosort2_default_parameters();
elseif options.KSversion==2.5
    job=UTkilosort2_default_parameters();
    job.sig        = 20;  % spatial smoothness constant for registration
    job.fshigh     = 300; % high-pass more aggresively
    job.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
    job.recompute_drift_map_after_shifting=1;
    job.kilosortVersion=2.5;
elseif options.KSversion==3
    job=UTkilosort3_default_parameters();
end

job.common_rejection_mode=options.common_rejection_mode;%'none';
job.remove_laser_artifact_sec=options.remove_laser_artifact_sec;%0;
job.interp_laser_sec=options.interp_laser_sec;%0;
job.rm_noise=options.rm_noise;
if isfield(options,'run_letters') && ~isempty(options.run_letters)
    job.runs_root=[BAPHYDATAROOT options.animal filesep options.site(1:end-1)];
else
    job.runs_root=[BAPHYDATAROOT options.animal filesep options.site(1:end-1)];
end
if isfield(options,'good_trials')
   if ~iscell(options.good_trials) || (length(options.good_trials)~=length(options.run_nums))
      error('options.good_trials must be cell array with length matching options.run_nums');
   end
   job.good_trials=options.good_trials;
end
job.skip_channels = eval(['[',options.skip_channels,']']);


if isfield(options,'run_letters') && ~isempty(options.run_letters)
    d=dir([job.runs_root filesep options.site(1:end-1) ,'*.m']);
else
    d=dir([job.runs_root filesep options.site ,'*.m']);
end
all_run_nums=cellfun(@(x)str2double(x(8:9)),{d.name});
all_site_letters=cellfun(@(x)x(7),{d.name});

dpsi=dir([job.runs_root filesep options.site ,'*/event_log.csv']);
psi_run_nums=cellfun(@(x)str2double(x((end-7):(end-6))),{dpsi.folder});
psi_site_letters=cellfun(@(x)x(end-8),{dpsi.folder});
for i = 1:length(psi_run_nums)
    [dpsi(i).name, dpsi(i).folder] = basename(dpsi(i).folder);
end

all_run_nums = [all_run_nums psi_run_nums];
all_site_letters = [all_site_letters psi_site_letters];
d = [d; dpsi];

if(ischar(options.run_nums))
    [~,si]=sort(all_run_nums);
    job.runs={d(si).name};
    options.run_nums=all_run_nums;
elseif isfield(options,'run_letters') && ~isempty(options.run_letters)
    for i=1:length(options.run_nums)
        ind=(options.run_nums(i)==all_run_nums) & (options.run_letters(i)==all_site_letters);
	        if(sum(ind)==0)
            error(['run ',num2str(options.run_nums(i)),' not found!'])
        end
        job.runs{i}=d(ind).name;
    end
else
    for i=1:length(options.run_nums)
        ind=options.run_nums(i)==all_run_nums;
        if(sum(ind)==0)
            error(['run ',num2str(options.run_nums(i)),' not found!'])
        end
        job.runs{i}=d(ind).name;
    end
end

if isfield(options,'whiten_split_runs')
    job.whiten_split_runs = options.whiten_split_runs;
end

%load first file to find datatype
dat=LoadMFile([job.runs_root filesep job.runs{1}]);

switch dat.globalparams.HWparams.DAQSystem
    case 'Open-Ephys'
        runinfo=rawgetinfo([job.runs_root filesep job.runs{1}],dat.globalparams);
        probe_id = getparm(options, 'probe','A');
        if strcmpi(probe_id,'')
            probe_id='A';
        end
        for i=1:length(runinfo.json_files)
            if strcmpi(runinfo.json_files(i).probe_id, probe_id)
                runinfo.json_file = runinfo.json_files(i).spikes;
                runinfo.json_file_lfp = runinfo.json_files(i).lfp;
                fprintf('probe_id=%s spikefile=%s\n', probe_id, runinfo.json_file);
            end
        end
        
        %Find and load channel map
        channel_map = getparm(options,'channel_map', '');
        if strcmpi(channel_map,'OE')
            [ports, probe_order_idx] = sort([runinfo.processors(1).NP_PROBE_info.port]);
            probe_number_idx = probe_order_idx(probe_id-'A'+1);
            NP_PROBE_INFO = runinfo.processors(1).NP_PROBE_info(probe_number_idx);
            fprintf('Probe%s: NP_PROBE .index=%d .port=%d .name=%s .serialno=%d\n', probe_id,...
                probe_number_idx, NP_PROBE_INFO.port, NP_PROBE_INFO.probe_name, NP_PROBE_INFO.probe_serial_number);
            % get channel map info from OE settings
            totalchannels=length(runinfo.processors(1).channels)
            cm = struct();
            cm.name='Neuropixels Phase3A';
            cm.probe_name=NP_PROBE_INFO.probe_name;
            cm.chanMap = [NP_PROBE_INFO.channels.number]+1;
            cm.chanMap0ind = [NP_PROBE_INFO.channels.number];
            cm.connected = zeros(1, totalchannels);
            cm.connected(runinfo.spike_channels)=1;
            cm.shankInd = ones(1, totalchannels);
            cm.xcoords = [NP_PROBE_INFO.channels.xpos];
            cm.ycoords = [NP_PROBE_INFO.channels.ypos];
            cm.bank = [NP_PROBE_INFO.channels.bank];
            job.chanMap = cm;
        elseif ~isempty(channel_map)
            %1st priority: If a channel_map was passed in options, use it
            fprintf('Using %s channel map as passed in options struct\n',channel_map)
            job.chanMap = [SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_' channel_map '.mat'];
        else
            %2nd priority: If a channel_map is filled out for this penetration in celldb, use it
            [channel_map, channel_map_path] = get_celldb_channel_map(options.site);
            if ~isempty(channel_map)
                fprintf('Using %s channel map specified in celldb\n',channel_map)
                job.chanMap = channel_map_path;
            elseif ~strcmp(runinfo.electrode,'unknown')
                %3rd priority: If a channel_map was given in the name of the channel mapping node in OEP, use it
                fprintf('Using %s channel map given by name of channel mapping node in OEP (saved in settings.xml)',runinfo.electrode)
                job.chanMap = [SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_' runinfo.electrode '.mat'];
                channel_map = runinfo.electrode;
            else
                if runinfo.spike_channels(1)==54 && length(runinfo.spike_channels)==64
                    channel_map = '64D_slot1';
                elseif runinfo.spike_channels(1)==11 && length(runinfo.spike_channels)==64
                    channel_map ='64D_slot1_bottom';    
                elseif runinfo.spike_channels(1)==118 && length(runinfo.spike_channels)==64
                    channel_map ='64D_slot2';
                elseif runinfo.spike_channels(1)==75 && length(runinfo.spike_channels)==64
                    channel_map ='64D_slot2_bottom';             
                elseif runinfo.spike_channels(1)==17 && length(runinfo.spike_channels)==128
                    channel_map = '128D_SepColsOffset';
                else
                    error(['Channel map not found in options struct, celldb penetration channel_map field, or by name of channel mapping node in OEP. ',...
                        'Could not be guessed based on channel map used in OEP signal chain.'])
                end
                warning(['Channel map not found in options struct, celldb penetration channel_map field, or by name of channel mapping node in OEP. ',...
                    'Guessed %s based on channel map used in OEP signal chain.'], channel_map)
                job.chanMap = [SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_' channel_map '.mat'];                
            end
        end
        if isstruct(job.chanMap)
            disp('extracted channel map from OE settings');
            ch=job.chanMap;
        elseif ~exist(job.chanMap,'file')
                error(['Channel map does not exist in %s. Create a channel map for this electrode.',...
            ' See /auto/users/lbhb/Code/multichannel_sorting/KiloSort2/configFiles/createChannelMapFile.m'])
            ch=load(job.chanMap,'chanMap');
        end
        job.probe_id = options.probe;
        job.Nchan = sum(ch.connected);
        if job.Nchan >= 384 %Neuropixels probes, use same # of clusters            
            job.Nfilt = job.Nchan;
        else % other probes use 1.5 times
            job.Nfilt = job.Nchan*1.5;
        end
        switch channel_map
            case {'neuropixPhase3A', 'neuropixPhase3Aall', 'OE'}
                job.keep_local = true; % flag to keep big local files for sorting and curating
            otherwise
                job.keep_local = false;
        end
        if isfield(options,'keep_local') && (job.keep_local ~= options.keep_local)
            fprintf('Overwriting keep_local from default of %d to %d as set in options.\n',job.keep_local,options.keep_local)
            job.keep_local = options.keep_local;
        end
        
    case 'MANTA'
        if dat.globalparams.NumberOfElectrodes==4
            fprintf('Four-channel MANTA recording detected, treating as a tetrode recording.\n')
            
            job.Nchan=4; %number of channels in the recording
            job.Nfilt= 32;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
            error('These param changes were for KiloSort1. See if anything should be changed for KiloSort2')
%             job.nNeighPC=4;
%             job.nNeigh=4;
%             job.whiteningRange=4;
%             job.spkTh=-4;
        else
            error(['This data looks like it was recorded with MANTA and has %d channels. ',...
                'Unknown channel mapping. You need to create a channel map for this electrode.',...
                ' See .../multichannel_sorting/KiloSort/configFiles/createChannelMapFile.m'],dat.globalparams.NumberOfElectrodes)
        end
    otherwise
        error('Unknown data format')
end
job.NchanTOT = job.Nchan;
job.datatype=dat.globalparams.HWparams.DAQSystem;

if options.Nfilt
   job.Nfilt=options.Nfilt;
end

job.name=options.site;
for i=1:length(options.run_nums)
    job.name=[job.name '_' num2str(options.run_nums(i))];
end

job.name=[job.name,'_KiloSort2'];

if(strcmp(job.common_rejection_mode,'none'))
    job.name=[job.name '_no_common_mode_rej'];
end

%% Drift correction options

%defaults for no drift correction
%job.find_drift_correction=false;
%job.fullMPMU_keep_spikes=false;
%job.keep_fproc=false; % keep in-progress Kilosort file
%job.driftCorrectionMode='None'; %no drift correction
%job.save_filtered_binary=false; 


%align similar clusters to peak as the same sample (helps avoid creation of
%two clusters with templated that look like offset copies of eachother)
%if 1
    %job.align_similar_clusters=true;
    %job.name=[job.name,'_alignSimilar'];
%end

% manually append ID to the same
if 0 % for testing which channel map
    job.name=[job.name,'_slot2'];
end

if getparm(options,'minfr_goodchannels',0)>0
    job.minfr_goodchannels=options.minfr_goodchannels;
    job.name=[job.name,sprintf('_minfr%.3f',options.minfr_goodchannels)];
else
    %job.name=[job.name,'_minfr_goodchannels_to0'];
    job.minfr_goodchannels=0;
end

if 0
    %mark bad channels
    job.chanMap='/auto/data/code/KiloSort/customs/chanMap_64D_slot1_SIK18_badchans_20by6to64.mat';
    job.chanMap='/auto/data/code/KiloSort/customs/chanMap_64D_slot1_bottom_badchans_center_column2.mat';
    ch=load(job.chanMap);
    job.Nchan=sum(ch.connected);
    job.NchanTOT=64;
    job.name=[job.name,'_rmCenterColumn'];
    
    if 0
       %to make a new channel map if other channels are bad: 
        ch=load('/auto/data/code/KiloSort/chanMap_64D_slot1');
        ch.connected(5)=0; %mark channel 5 bad. Change to whatever is bad
        new_name='/auto/data/code/KiloSort/chanMap_64D_slot1_badchans_blah_blah'; % change to say which channels are marked bad
        save(new_name,'-Struct','ch') % save it
        %Then use the new map in job.chanMap above
    end
end
if job.kilosortVersion==3
    job.name=strrep(job.name,'KiloSort2','KiloSort3');
    matlab_str='UTkilosort3_run_job';
elseif job.kilosortVersion==2.5
    job.name=strrep(job.name,'KiloSort2','KiloSort2pt5');
    matlab_str='UTkilosort2pt5_run_job';
else
    matlab_str='UTkilosort2_run_job';
end

if 0 && options.KSversion>=2.5 % this breakd KS3
    job.name=[job.name,'_rigid'];
    job.nblocks=1; % why is this rigid??
end

if 0 && options.KSversion>=2.5
    job.name=[job.name,'_batchsize256'];
    job.NT = 256*1024+ job.ntbuff; % batch size , must be multiple of 32 + ntbuff. Reduce size if running out of memory
end

if 0 && options.KSversion>=2.5
    job.name=[job.name,'_medfilt3'];
    job.drift_estimate_median_filter_pts=3;
end

if isfield(options,'Th') & (sum(abs(options.Th-job.Th))>0)
    job.name=[job.name,sprintf('_Th%d_%d',options.Th(1),options.Th(2))];
    job.Th=options.Th;
end

if isfield(options,'AUCsplit') && (options.AUCsplit ~= job.AUCsplit)
    %If an AUCsplit was passed and it's not equal to the default, change it and append it to the job name
    job.name=[job.name,'_AUCsplit',sprintf('%.0f',options.AUCsplit*100)];
    job.AUCsplit=options.AUCsplit;
end

if isfield(options,'probe') && ~strcmpi(options.probe,'')
    %if probe specified (simultaneous neuropix) add to name
    job.name=[job.name,'_Probe',options.probe];
end
if isfield(options,'rm_noise') && (options.rm_noise==1)
    %if probe specified (simultaneous neuropix) add to name
    job.name=[job.name,'_rmnoise'];
    job.rm_noise=options.rm_noise;
end

job.results_path=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'results'];
job.results_path_temp=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'big'];
job.fbinary=[job.results_path_temp filesep 'binary.dat'];% will be created for 'openEphys'
job.fproc=[job.results_path_temp filesep 'temp_wh.dat'];% residual from RAM of preprocessed data
job.status=0; %0: not started, 1:sorted, 2: manually analyzed, 3: completed
disp(job)
fn=[KSpaths.jobs_path filesep 'in_progress' filesep job.name '.mat'];


if(exist(fn,'file'))
    if isfield(options,'allow_overwrite') && options.allow_overwrite
        fprintf('Deleted existing job file, %s. Job file (%s) existed, and options.allow_overwrite was True.\n',fn)
        delete(fn)
    else
        error(['This job already exists! Filepath: ',fn,' set options.delete_and_overwrite_if_job_exists to True to delete it'])
%         f = questdlg(['This job already exists! Filepath: ',fn,...
%             ' Delete and overwrite or cancel?'],'File Exists','Overwrite','Cancel','Overwrite');
%         if strcmp(f,'Overwrite')
%             delete(fn)
%         else
%             fprintf('\nCancelled')
%             return
%         end
    end
end

UTmkdir(fn);
save(fn,'-Struct','job');
if ispc
    error('Need to open permissions...Figure out how to set permissions for all users from a Windows machine')
else
    [w,s]=unix(['chmod 777 ',fn]);if w, error(s), end
end
fprintf(['Created job: ',fn,'\n'])

if job.kilosortVersion==3
    %UTkilosort3_run_job(fn);
    UTkilosort2_run_job(fn);
elseif job.kilosortVersion==2.5
    UTkilosort2_run_job(fn);
else
    UTkilosort2_run_job(fn);
end
