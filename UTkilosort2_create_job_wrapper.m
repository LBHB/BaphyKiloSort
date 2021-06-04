function UTkilosort2_create_job_wrapper(options)

%script for creating KiloSort jobs and adding them to the jobs folder
global BAPHYDATAROOT SERVER_PATH
if(isempty(BAPHYDATAROOT))
    baphy_set_path
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
        % override electrode name, if specified in options
        electrode_name = getparm(options,'channel_map', runinfo.electrode);
        switch electrode_name
            case {'64D_slot1','64D_slot2','64M_slot1','64M_slot2', '64D_slot1_bottom', '64D_slot2_bottom'}
                job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_' electrode_name '.mat'];
                job.Nchan=64;
                job.Nfilt= 96;                
            case {'128D'}
                [s, UCLA_to_OEP]=probe_128D();
            case {'128P_bottom'}
                job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_' electrode_name '.mat'];
                job.Nchan=128;
                job.Nfilt=192;
            case 'unknown'
                if runinfo.spike_channels(1)==54 && length(runinfo.spike_channels)==64
                    %channel map for UCLA 64D in Intan headstage slot 1
                    job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_64D_slot1.mat'];
                    job.Nchan=64;
                    job.Nfilt= 96;
                    ename='64D_slot1';
                elseif runinfo.spike_channels(1)==11 && length(runinfo.spike_channels)==64
                    job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_64D_slot1_bottom.mat'];
                    job.Nchan=64;
                    job.Nfilt= 96;
                    ename='64D_slot1_bottom';
                elseif runinfo.spike_channels(1)==118 && length(runinfo.spike_channels)==64
                    %channel map for UCLA 64D in Intan headstage slot 2
                    job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_64D_slot2.mat'];
                    job.Nchan=64;
                    job.Nfilt= 96;
                    ename='64D_slot2';
                elseif runinfo.spike_channels(1)==17 && length(runinfo.spike_channels)==128
                    %channel map for UCLA 128D w/ Intan headstage
                    job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_128D_SepColsOffset.mat'];
                    job.Nchan=128;
                    job.Nfilt= 192;
                    ename='128D';
                else
                    error(['Electrode not specified in Settings.xml, and can''t be inferred by chappel mapping (not 6D, 64B, or 128D).'...
                        ' Rename channel map module in Open Ephys to save electrode name. Then create a channel map for this electrode and add it to this code. ',...
                        ' See .../multichannel_sorting/KiloSort/configFiles/createChannelMapFile.m'])
                end
                warning(['Electrode was not specified in Settings.xml. Rename channel map module in Open Ephys to save electrode name.'...
                    ' Guessed ' ename 'based on channel map.'])
            otherwise
                error(['Open-Ephys settings.xml said that the electrode was: "' electrode_name '". No Kilosort channel map file exists for this electrode.'...
                    'Create a channel map for this electrode and add it to this code. ',...
                        ' See .../multichannel_sorting/KiloSort/configFiles/createChannelMapFile.m'])
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
    job.name=[job.name,'_minfr_goodchannels_to0'];
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

job.name=[job.name ,'_iNeigh_LASv2'];

job.results_path=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'results'];
job.results_path_temp=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'big'];
job.fbinary=[job.results_path_temp filesep 'binary.dat'];% will be created for 'openEphys'
job.fproc=[job.results_path_temp filesep 'temp_wh.dat'];% residual from RAM of preprocessed data
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

if(1)
    rundataid=0;
    if job.kilosortVersion==3
        progname='matlabbgLAS queuerun';
        parmstring=['UTkilosort3_run_job(''',fn,''');'];
    elseif job.kilosortVersion==2.5
        progname='matlabbgLAS queuerun';
        parmstring=['UTkilosort2pt5_run_job(''',fn,''');'];
    else
        progname='matlabbg queuerun';
        parmstring=['UTkilosort2_run_job(''',fn,''');'];
    end
    allowqueuemaster=1;
    note=['Kilosort job: ',job.name];
    user=getenv('USER');
    user = 'nems';
    r=dbGetConParms();
    rtemp=r;
    rtemp.DB_SERVER='hyrax.ohsu.edu:3306';
    
    % temporarily connect to NEMS db
    dbSetConParms(rtemp);
    
    if ~exist('dbaddqueue','file')
        narf_set_path
    end
    queueidx=dbaddqueue(rundataid,progname,parmstring,allowqueuemaster,note,user,job.GPU);
    save(fn,'-append','queueidx');
    job.queueidx=queueidx;
    
    % restore old db connection
    dbSetConParms(r);
    
    fprintf('\n')
    disp(['<a href="http://hyrax.ohsu.edu/celldb/queuemonitor.php?user=%25&complete=-2&machinename=%25&notemask=',job.name,'">Check job status</a>'])
    disp(['<a href="http://hyrax.ohsu.edu/celldb/queuemonitor.php?user=%&complete=-1&machinename=hyena&notemask=">See what''s running on Hyena</a>'])
end    