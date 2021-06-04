function ops=UTkilosort2_default_parameters

ops.kilosortVersion = 2;
ops.verbose             = 1; % whether to print command line progress		
		
ops.datatype            = 'openEphys';  % binary ('dat', 'bin') or 'openEphys'		
%ops.fbinary             = fullfile(fpath, 'binary.dat'); % will be created for 'openEphys'		
%ops.fproc               = fullfile(fpath, 'temp_wh.dat'); % residual from RAM of preprocessed data		
%ops.root                = {'/tmp/evpread/Boleto/BOL007/raw/BOL007d01/BOL007d01_p_BNB_2016-10-17_12-17-00','/tmp/evpread/Boleto/BOL007/raw/BOL007d02/BOL007d02_p_FTC_2016-10-17_12-18-49/'}; % 'openEphys' only: where raw files are		
% define the channel map as a filename (string) or simply an array		
%ops.chanMap             = fullfile('/auto/users/luke/Code/KiloSort_baphy/chanMap.mat'); % make this file using createChannelMapFile.m		
ops.common_rejection_mode = 'none';  % 'mean' 'median'
ops.raw_scale_factor=1;


%% Kilosort2 ops
ops.trange = [0 Inf]; % time range to sort

ops.fig=0; % show figures during template matching

% sample rate
ops.fs = 30000;  

% frequency for high pass filtering (150)
ops.fshigh = 150;   

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 10;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 
%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

