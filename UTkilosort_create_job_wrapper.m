function UTkilosort_create_job_wrapper(options)

%script for creating KiloSort jobs and adding them to the jobs folder
global BAPHYDATAROOT SERVER_PATH
if(isempty(BAPHYDATAROOT))
    baphy_set_path
end

options.animal=getparm(options,'animal','Tartufo');
options.site=getparm(options,'site','TAR017b');
options.run_nums=getparm(options,'run_nums',[10]);
options.Nfilt=getparm(options,'Nfilt',0);
options.common_rejection_mode=getparm(options,'common_rejection_mode','median');


KSpaths=UTkilosort_paths();


% All this code is to fill in job.runs_root and job.runs
job=UTkilosort_default_parameters();
job.kilosortVersion = 1;
job.common_rejection_mode=options.common_rejection_mode;%'none';
job.runs_root=[BAPHYDATAROOT options.animal filesep options.site(1:end-1)];
job.GPU=1;
job.showfigures=0;
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
        if runinfo.spike_channels(1)==54 && length(runinfo.spike_channels)==64
            %channel map for UCLA 64D in Intan headstage slot 1
            job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_64D_slot1.mat'];
            job.Nchan=64;
            job.Nfilt= 96;
        elseif runinfo.spike_channels(1)==118 && length(runinfo.spike_channels)==64
            %channel map for UCLA 64D in Intan headstage slot 2
            job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_64D_slot2.mat'];
            job.Nchan=64;
            job.Nfilt= 96;
        elseif runinfo.spike_channels(1)==17 && length(runinfo.spike_channels)==128
            %channel map for UCLA 128D w/ Intan headstage
            job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_128D_SepColsOffset.mat'];
            job.Nchan=128;
            job.Nfilt= 192;
        else
            error(['Unknown channel mapping. You need to create a channel map for this electrode.',...
                ' See .../multichannel_sorting/KiloSort/configFiles/createChannelMapFile.m'])
        end
    case 'MANTA'
        if dat.globalparams.NumberOfElectrodes==4
            fprintf('Four-channel MANTA recording detected, treating as a tetrode recording.\n')
            job.Nchan=4; %number of channels in the recording
            job.Nfilt= 32;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
            job.nNeighPC=4;
            job.nNeigh=4;
            job.whiteningRange=4;
            job.spkTh=-4;
        else
            error(['This data looks like it was recorded with MANTA and has %d channels. ',...
                'Unknown channel mapping. You need to create a channel map for this electrode.',...
                ' See .../multichannel_sorting/KiloSort/configFiles/createChannelMapFile.m'],dat.globalparams.NumberOfElectrodes)
        end
    otherwise
        error('Unknown data format')
end
job.datatype=dat.globalparams.HWparams.DAQSystem;

%job.runs_root and job.runs could also be filled manually if desired:
%job.runs_root=[BAPHYDATAROOT 'Boleto' filesep 'BOL007'];
%job.runs={'BOL007e01_p_BNB','BOL007e02_p_FTC','BOL007e03_p_TOR','BOL007e04_p_FTC'};

%These are now defined above. Overwrite here if you want.
%job.Nchan=64; %number of channels in the recording
%job.Nfilt= 96;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
%job.chanMap=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'chanMap_128D.mat']; %channel map for UCLA 128D w/ Intan headstage
if options.Nfilt
   job.Nfilt=options.Nfilt;
end

job.name=options.site;
for i=1:length(options.run_nums)
    job.name=[job.name '_' num2str(options.run_nums(i))];
end
job.name=[job.name,'_',num2str(job.Nfilt),'clusts'];
if(strcmp(job.common_rejection_mode,'none'))
    job.name=[job.name '_no_common_mode_rej'];
end

%% Drift correction options

%defaults for no drift correction
job.find_drift_correction=false;
job.fullMPMU_keep_spikes=false;
job.keep_fproc=false; % keep in-progress Kilosort file
job.driftCorrectionMode='None'; %no drift correction
job.save_filtered_binary=false; 
 
%Turn on various drift-correction things (in order of most useful)
if 0
    % save filtered binary if it will be used for drift correction. Can also be created later
    job.save_filtered_binary=true;
elseif 0
    %Apply drift correction file
    job.driftCorrectionMode='BeforeFiltering';
    %job.driftCorrectionMode='AfterFiltering';
    job.driftCorrectionInterpMode='1d';
    %job.driftCorrectionInterpMode='2d';
    %try changing parameters, doesn't seem to help much
    % job.Th               = [3 8 8];
    % job.lam              = [7 7 7];
    %job.mergeT           = .2;
    %job.splitT           = .2;
    
    %define which drift correction file to use
    %job.driftCorrectionFile=['/auto/users/luke/Projects/MultiChannel/Drift/images/',job.name,'/size 6 step 0.6.mat'];
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre196b_4_5_6_7_96clusts_keep_fproc_and_spikes/size 12 step 0.6 windowed allspikes bounded clustered weight_1ov4_lt50.mat';
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre196b_4_5_6_7_96clusts_keep_fproc_and_spikes/size 12 step 0.6 windowed allspikes bounded clustered weight_1ov4_lt50 shift_p30sec.mat';
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre196b_9_10_11_12_96clusts/size 60 step 1 amp_60 windowed allspikes bounded clustered weight_1ov4_lt50.mat';
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre197c_1_2_3_192clusts/size 20 step 1 amp_60 windowed allspikes bounded clustered weight_1ov4_lt50 nomspc.mat';
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre197c_1_2_3_4_5_6_7_8_9_10_11_12_192clusts_alignSimilar/size 60 step 1 amp_60 windowed allspikes bounded clustered weight_1ov4_lt50 nomspc.mat';
    job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre197c_1_2_3_4_5_6_7_8_9_10_11_12_192clusts_alignSimilar/size 30 step 1 KIDs_112  130 windowed allspikes bounded clustered weight_1ov4_lt50 nomspc.mat';
    %job.driftCorrectionFile='/auto/users/luke/Projects/MultiChannel/Drift/images/fre196b_6_96clusts_Alpaca/size 6 step 6.000000e-01 2x.mat';
    
    %change job name 
    job.name=[job.name,'_corrected_by_clustered_before_filtering_30_KIDs112_113'];
elseif 0
    %Create and apply drift corrction before filtering
    job.find_drift_correction=true;
    job.driftCorrectionMode='BeforeFiltering';
    job.name=[job.name,'_corrected_before_filtering'];
    job.driftCorrectionFile=['/auto/users/luke/Projects/MultiChannel/Drift/images/',job.name,'/size 6 step 0.6.mat'];
elseif 0
    % For testing cortex-lab (Marius') drift correction (didn't work)
    job.driftCorrectionMode='cortex-lab';
    job.name=[job.name,'_dc_cortex-lab_V2'];
elseif 0
    % Keep all spikes to be used to calculate drift correction later (too big for big jobs)
    job.fullMPMU_keep_spikes=true;
    job.name=[job.name,'_keep_spikes'];
    job.find_drift_correction=false;
end

%align similar clusters to peak as the same sample (helps avoid creation of
%two clusters with templated that look like offset copies of eachother)
if 1
    job.align_similar_clusters=true;
    job.name=[job.name,'_alignSimilar'];
end

% manually append ID to the same
if 0 % for testing which channel map
    job.name=[job.name,'_slot2'];
end

job.results_path=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'results'];
job.results_path_temp=[job.runs_root filesep 'tmp' filesep 'KiloSort' filesep  job.name filesep 'big'];
job.fbinary=[job.results_path_temp filesep 'binary.dat'];% will be created for 'openEphys'
job.fproc=[job.results_path_temp filesep 'temp_wh.dat'];% residual from RAM of preprocessed data
job.status=0; %0: not started, 1:sorted, 2: manually analyzed, 3: completed

fn=[KSpaths.jobs_path filesep 'in_progress' filesep job.name '.mat'];
if(exist(fn,'file'))
    delete(fn)
end

if(exist(fn,'file'))
    error(['This job already exists! Filepath: ',fn])
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
    progname='matlabbg queuerun';
    parmstring=['UTkilosort_run_job(''',fn,''');'];
    allowqueuemaster=1;
    note=['Kilosort job: ',job.name];
    user=getenv('USER');
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