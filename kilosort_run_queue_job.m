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
options.run_nums=getparm(options,'run_nums',[10]); %can be 'all'
options.Nfilt=getparm(options,'Nfilt',0);
options.common_rejection_mode=getparm(options,'common_rejection_mode','median');
options.KSversion=getparm(options,'KSversion',2);
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
    job.recompute_drift_map_after_shifting=1;
    job.kilosortVersion=2.5;
elseif options.KSversion==3
    job=UTkilosort3_default_parameters();
end

job.common_rejection_mode=options.common_rejection_mode;%'none';
job.remove_laser_artifact_sec=options.remove_laser_artifact_sec;%0;
job.interp_laser_sec=options.interp_laser_sec;%0;
job.runs_root=[BAPHYDATAROOT options.animal filesep options.site(1:end-1)];
if isfield(options,'good_trials')
   if ~iscell(options.good_trials) || (length(options.good_trials)~=length(options.run_nums))
      error('options.good_trials must be cell array with length matching options.run_nums');
   end
   job.good_trials=options.good_trials;
end

d=dir([job.runs_root filesep options.site ,'*.m']);
all_run_nums=cellfun(@(x)str2double(x(8:9)),{d.name});
if(ischar(options.run_nums))
    [~,si]=sort(all_run_nums);
    job.runs={d(si).name};
    options.run_nums=all_run_nums;
else
    for i=1:length(options.run_nums)
        ind=options.run_nums(i)==all_run_nums;
        if(sum(ind)==0)
            error(['run ',num2str(options.run_nums(i)),' not found!'])
        end
        job.runs{i}=d(ind).name;
    end
end

%load first file to find datatype
dat=LoadMFile([job.runs_root filesep job.runs{1}]);
switch dat.globalparams.HWparams.DAQSystem
    case 'Open-Ephys'
        runinfo=rawgetinfo([job.runs_root filesep job.runs{1}],dat.globalparams);
        
        %Find and load channel map
        channel_map = getparm(options,'channel_map', '');
        if ~isempty(channel_map)
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
        if ~exist(job.chanMap,'file')
                error(['Channel map does not exist in %s. Create a channel map for this electrode.',...
            ' See /auto/users/lbhb/Code/multichannel_sorting/KiloSort2/configFiles/createChannelMapFile.m'])
        end
        ch=load(job.chanMap,'chanMap');
        job.Nchan = length(ch.chanMap);
        if job.Nchan >= 384 %Neuropixels probes, use same # of clusters            
            job.Nfilt = job.Nchan;
        else % other probes use 1.5 times
            job.Nfilt = job.Nchan*1.5;
        end
        switch channel_map
            case 'neuropixPhase3A'
                job.keep_local = true; % flag to keep big local files for sorting and curating
            otherwise
                job.keep_local = false;
        end
        if isfield(options,'keep_local')
            fprintf('Ovverwriting keep_local from default of %d to %d as set in options.\n',job.keep_local,options.keep_local)
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

if 1
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

if 1 && options.KSversion>=2.5
    job.name=[job.name,'_rigid'];
    job.nblocks=1;
end

if 1 && options.KSversion>=2.5
    job.name=[job.name,'_batchsize256'];
    job.NT = 256*1024+ job.ntbuff;
end

if 1 && options.KSversion>=2.5
    job.name=[job.name,'_medfilt3'];
    job.drift_estimate_median_filter_pts=3;
end

if 0 
    job.name=[job.name,'_Th10_2'];
    job.Th=[10 2];
end

job.results_path=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'results'];
job.results_path_temp=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'big'];
job.fbinary=[job.results_path_temp filesep 'binary.dat'];% will be created for 'openEphys'
job.status=0; %0: not started, 1:sorted, 2: manually analyzed, 3: completed

fn=[KSpaths.jobs_path filesep 'in_progress' filesep job.name '.mat'];


if(exist(fn,'file'))
    f = questdlg(['This job already exists! Filepath: ',fn,...
        ' Delete and overwrite or cancel?'],'File Exists','Overwrite','Cancel','Overwrite');
    if strcmp(f,'Overwrite')
        delete(fn)
    else
        fprintf('\nCancelled')
        return
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
    UTkilosort3_run_job(fn);
elseif job.kilosortVersion==2.5
    UTkilosort2pt5_run_job(fn);
else
    UTkilosort2_run_job(fn);
end
