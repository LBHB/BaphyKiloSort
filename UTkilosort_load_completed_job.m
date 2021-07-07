function savfiles=UTkilosort_load_completed_job(job_file,options)
%UTkilosort_load_completed_job(job_file)
%UTkilosort_load_completed_job(job_file,options)
% converts templates and spike assignments created by kilosort (and phy
% if present) to *.spk.mat file readable by Baphy.

% Inputs:
% job_file: path to file containing parameters for this job
% options: structure of options:
    % load_as_temp: if true, stores results in temporary location instead of on server
    % force_compute_quality_measures: if true, computes quality measures (these
    % will be computed by default if load_ as_temp is true, but not if it is
    % false).

% Outputs:
% savefiles: paths to created .spk.mat file
%
% Comments on system for saving to celldb (SVD 2018-07-02)
%
% calls dblink_phy to save meta data associated with each spk.mat structure to 
%  the corresponding sCellFile entries in celldb
% if append_units=True, data will be assigned to NEW units, starting with
% the next unit number following the current highest value. IMPORTANT: if 
% you re-sort and save part of the recording and have used append_units,
% you should probably delete all the relevant celldb data associated with 
% the current site and then resave everything. 
%
% to purge, run these mysql commands (which <siteid> is the relevant siteid):
%   delete from gSingleRaw where cellid like "<siteid>%";
%   delete from gSingleCell where cellid like "<siteid>%";
%   delete from sCellFile where cellid like "<siteid>%";
%
fprintf(['\nSaving ',job_file,'\n'])
if ~exist('readNPY.m','file')
    error('readNPY not found. Is MULTICHANNEL_SORTING_PATH set so that the npy-matlab package is available?')
end

options.use_automerge=getparm(options,'use_automerge',0); %Depreciated (sued by Kilosort 1 but not > 2)
options.load_as_temp=getparm(options,'load_as_temp',0);
options.force_compute_quality_measures=getparm(options,'force_compute_quality_measures',0);
options.append_units=getparm(options,'append_units',0);
options.delete_existing_spkfile=getparm(options,'delete_existing_spkfile',options.load_as_temp);
options.purge_rawids=getparm(options,'purge_rawids',[]);
options.remove_Kids=getparm(options,'remove_Kids',[]); %Kid means KiloSort id
options.keep_Kids=getparm(options,'keep_Kids',[]);

%Options currently not set by Kilosort_browser
options.unit_numbering_offset=getparm(options,'unit_numbering_offset',0); %use this to make units 1:3 be units 11:13 (for example if you want to make the unique from other runs)
options.merge_all_units=getparm(options,'merge_all_units',0); %Merges all units into a single unit, assigned to the channel of the lowest Kid number
options.force_best_channels_to_one=getparm(options,'force_best_channels_to_one',0); %Sets all units to be on channel 1

if ~isempty(options.keep_Kids) && ~isempty(options.remove_Kids)
    error('Either define keep_Kids to tell me which ids to keep, or remove_Kids to tel me which to remove')
end
    
if(options.use_automerge)
    suffix='_after_automerge';
else
    suffix='';
end
job=load(job_file);
for i=1:length(job.runs)
    LoadMFile([job.runs_root filesep job.runs{i}]);
    globalparams_{i}=globalparams;
    exptparams_{i}=exptparams;
    exptevents_{i}=exptevents;
end

if(options.use_automerge)
    rez=load([job.results_path_temp,filesep,'rez.mat'],'Wraw');
    rez_am=load([job.results_path_temp,filesep,'rez',suffix,'.mat'],'st3');
    rez.st3=rez_am.st3;
    clear rez_am
elseif isfield(job,'kilosortVersion') && job.kilosortVersion>=2
    rez=load([job.results_path_temp,filesep,'rez',suffix,'.mat'],'st3','W','U','mu', 'ops');
    % Kilosort2 doesn't create Wraw (templates in original space) at the 
    % end of the algorithm so do it here.
    rez.Wraw = [];
    % Kilosort2 may remove channels that don't mean min FR criterion.
    % check if this is the case
    ch_map = load(job.chanMap);
    chanMap = ch_map.chanMap(ch_map.connected);
    kept_chans = rez.ops.chanMap;
    if size(kept_chans, 1) < size(chanMap, 1)
        if ~isempty(strfind(job.chanMap, 'slot2'))
            rm_chans = setdiff(chanMap-64, kept_chans);
        else
            rm_chans = setdiff(chanMap, kept_chans);
        end
        % then, remap back to electrode nums
        removed_inds = zeros(size(rm_chans));
        for i=1:length(rm_chans)
            removed_inds(i)=find(rm_chans(i)==(chanMap'-64));
        end
        kept_inds = 1:length(chanMap);
        kept_inds(ismember(kept_inds, removed_inds)) = [];
    end
    
    for n = 1:size(rez.U,2)
        % temporarily use U rather Urot until I have a chance to test it
        rez.Wraw(:,:,n) = rez.mu(n) * squeeze(rez.U(:,n,:)) * squeeze(rez.W(:,n,:))';
    end
    
    %Spike are not sorted by time like they were in KS1, sort them here
    [~,si]=sort(rez.st3(:,1));
    rez.st3=rez.st3(si,:);
    
else
    rez=load([job.results_path_temp,filesep,'rez',suffix,'.mat'],'st3','Wraw');
end
clusts=readNPY([job.results_path,suffix,filesep,'spike_clusters.npy']);
template_ids=readNPY([job.results_path,suffix,filesep,'spike_templates.npy']);
%templates_unwitened=readNPY([job.results_path,'templates_unw.npy']); %(64x82x64) % unwhitened templates
%templates=readNPY([job.results_path,'templates.npy']); %(64x82x64) % templates

%a=readNPY([job.results_path,'pc_features.npy']); %PCs (N*3*12) why 12?:
%for 12 different most-informaive channels, channels used stored in pc_feature_ind
%b=readNPY([job.results_path,'pc_feature_ind.npy']);

cluster_group_file=[job.results_path,suffix,filesep,'cluster_groups.tsv'];
if(~exist(cluster_group_file,'file'))
    cluster_group_file=[job.results_path,suffix,filesep,'cluster_group.tsv'];
    if(~exist(cluster_group_file,'file'))
        error('Cluster group file (%s) does not exist and is needed to load completed jobs. Did you save in phy?',cluster_group_file)
    end
end
fid= fopen(cluster_group_file);
C = textscan(fid,'%d %s',-1,'Headerlines',1,'Delimiter','\t');
fclose(fid);
ids=C{1};
group=C{2};
if(options.merge_all_units)
    ids=min(ids);
    group={'mua'};
    clusts(:)=min(clusts);
end

% Keep only group == 'good' | group == 'mua' and update ids
su_mu_mask = cellfun(@(x) (strcmp(x,'good')|strcmp(x, 'mua')), group);
group = group(su_mu_mask);
ids = ids(su_mu_mask);

sus=find(strcmp(group,'good'));
mus=find(strcmp(group,'mua'));
%sus=find(strcmp(group,'unsorted'));

all_units=sort([sus;mus]); %units to save
%all_units=sort([sus]); %units to save

if ~isempty(options.remove_Kids)
    mi = ismember(ids(all_units),options.remove_Kids);
    if sum(mi) < length(options.remove_Kids)
        s1=sprintf('%d ', options.remove_Kids);s1(end)=[];
        s2=sprintf('%d ', ids(all_units)); s2(end)=[];
        a=questdlg(['Warning! Requested ids to remove (',s1,') not all found in id list (',s2,')'],'Dont''t save ids not found','Continue','Cancel','Cancel');
        if strcmp(a,'Cancel')
            savfiles=[];
            fprintf('\n Cancelled')
            return
        end
        end
    all_units(ismember(ids(all_units),options.remove_Kids))=[];
end


if(isempty(all_units))
    error('No units found, not saving')
end
best_chan_file=[job.results_path,suffix,filesep,'best_channels.npy'];
if(exist(best_chan_file,'file'))
    best_channels_phy=readNPY(best_chan_file)+1;
else
    best_channels_phy=[];
end
if(options.merge_all_units)
    best_channels_phy=1;
end
snrs_file=[job.results_path,suffix,filesep,'snrs.npy'];
if(exist(snrs_file,'file'))
    snrs=readNPY(snrs_file);
else
    snrs=[];
end

mean_waveforms_file=[job.results_path,suffix,filesep,'mean_waveforms.npy'];
if(exist(mean_waveforms_file,'file'))
    mean_waveforms_phy=readNPY(mean_waveforms_file);
else
    mean_waveforms_phy=[];
end

% Reading in wave stats and waveform classifications
% CRH 2018-06-15
wft_struct = struct();
all_files = dir(fullfile([job.results_path,suffix], 'wft_*'));
for i = 1:length(all_files)
    n = all_files(i).name;
    value = readNPY([job.results_path,suffix,filesep,n]);
    fparts = strsplit(n, '.');
    field = fparts{1};
    wft_struct.(field) = value;
end
if ~(isempty(fieldnames(wft_struct))) && isfield(wft_struct, 'wft_mwf')
    wft_struct = rmfield(wft_struct,'wft_mwf');  % no need to save mean waveform
end

% will use later to define cellids
fparts = strsplit(job.results_path, filesep);
site = fparts(9);
site = site{1};
site = site(1:7);

merge_results_file=[job.results_path,suffix,filesep,'cluster_names.ts'];
if ~options.append_units && (exist(merge_results_file,'file'))
    fprintf('Merge results file exists, using: %s\n',merge_results_file)
    fid= fopen(merge_results_file);
    C = textscan(fid,'%d %d %d %d',-1,'Headerlines',1,'Delimiter','\t');
    fclose(fid);
    ids_phy=C{1};
    unit_types_phy=C{2};
    best_channels_phy_merge=C{3}+1;
    unit_numbers_phy=C{4};
    
    %remove noise clusters
    rmi=find(unit_types_phy==3);
    ids_phy(rmi)=[];
    unit_types_phy(rmi)=[];
    best_channels_phy_merge(rmi)=[];
    unit_numbers_phy(rmi)=[];
    
    if(~isempty(snrs))
        [~,ind]=ismember(ids_phy,ids);
        snrs=snrs(:,ind);
    end
    
    % copying what's being done for snrs CRH
    if(~isempty(wft_struct))
        [~, ind]=ismember(ids_phy,ids);
        fn = fieldnames(wft_struct);
        for i = 1:length(fn)
            key = char(fn(i));
            value = wft_struct.(key);
            if strcmp(fn(i), 'wft_mwfs')
                wft_struct.(key) = value(:, ind);
            else
                wft_struct.(key) = value(ind);
            end
        end
    end
    
    ids=ids_phy;
    all_units=1:length(ids);
else
    best_channels_phy_merge=[];
end


%st3: times (in samples), templates, amplitudes
spiketimes=rez.st3(:,1);

%spikeTemplates = uint32(rez.st3(:,2));
spike_assignments=clusts;

best_channel=zeros(length(all_units),1,'int8');
was_split=false(length(all_units),1);
was_merged=false(length(all_units),1);
for i=1:length(all_units)
    [templates(:,:,i),best_channel(i),was_split(i),was_merged(i)] = UTkilosort_get_template(clusts,template_ids,ids(all_units(i)),rez);
    if(i==1)
        templates(:,:,2:length(all_units))=NaN;
    end
end

if(~isempty(best_channels_phy_merge))
    best_channel=best_channels_phy_merge;
elseif(~isempty(best_channels_phy))
    %if phy saved best channels, use these instead of what's computed in UTkilosort_get_template
    best_channel=best_channels_phy(all_units);
end

if(options.force_best_channels_to_one)
    best_channel(:)=1;
end

spiketimes_by_trial=nan(3,length(spiketimes));
switch job.runinfo(1).evpv
    case 5
        recording_start=[0 cumsum(job.nSamplesBlocks(1:end-1))]+1;
        trial_onsets_= job.trial_onsets_;
        runs_per_trial=cellfun('length',trial_onsets_);
        for i=1:length(job.runs)
            trial_onsets_{i}=trial_onsets_{i}+sum(job.nSamplesBlocks(1:i-1));
        end
    case 6 %6 means Open-Ophys format
        if isfield(job,'runs_per_trial')
            trial_onsets_=job.trial_onsets_;
            runs_per_trial = job.runs_per_trial;
        else
            str='This job file looks old (it doesn''t have the field ''runs_per_trial''). If it was recorded with open ephys 0.4.5 there are probably timing errors. You can fix them by running UTkilosort2_fix_trial_onsets. If it''s older, it might not have these errors so you could probably safely skip this step.';
            resp = questdlg(str,'Old File, fix trial onsets?','run UTkilosort_fix_trial_onsets','use trial onsets in all_channels.events','run UTkilosort_fix_trial_onsets');
            switch resp
                case 'run UTkilosort_fix_trial_onsets'
                    job=UTkilosort2_fix_trial_onsets(job_file);
                    trial_onsets_=job.trial_onsets_;
                    runs_per_trial = job.runs_per_trial;                    
                case 'use trial onsets in all_channels.events'
                    for i=1:length(job.runs)
                        not_filesep=setdiff({'/','\'},filesep);not_filesep=not_filesep{1};
                        globalparams_{i}.rawfilename=strrep(globalparams_{i}.rawfilename,not_filesep,filesep);
                        
                        [~, EVinfo(i),EVtimestamps{i}] = load_open_ephys_data_faster([globalparams_{i}.rawfilename,filesep,'all_channels.events']);
                        recording_start(i)=EVtimestamps{i}(1);
                        trial_onsets_{i}=(EVtimestamps{i}(EVinfo(i).eventId==1&EVinfo(i).eventType==3)-recording_start(i)+1)'+sum(job.nSamplesBlocks(1:i-1)); % in samples re recording start
                        runs_per_trial(i)=length(trial_onsets_{i});
                        sampleRate(i)=EVinfo(i).header.sampleRate;
                    end
                otherwise
                    savfiles=[];
                    return
            end
                    
        end
end
trial_onsets=[trial_onsets_{:}];

%add recording length to trial onsets(2:end)

run_starts=[1 cumsum(runs_per_trial)+1];
trial_onsets(end+1)=inf; %add inf to store any spikes after last trial at end of last trial
spike_start_idx_by_run=nan(size(runs_per_trial));
for trialidx=1:length(trial_onsets)-1
    idx=find(spiketimes>=trial_onsets(trialidx) & spiketimes<trial_onsets(trialidx+1));
    this_run=find(run_starts<=trialidx,1,'last');
    if(isnan(spike_start_idx_by_run(this_run)) && ~isempty(idx))
        spike_start_idx_by_run(this_run)=idx(1);
    end
    trialidx_per_run=trialidx-sum(runs_per_trial(1:this_run-1));
    spiketimes_by_trial(1,idx)=trialidx_per_run;
    spiketimes_by_trial(2,idx)=spiketimes(idx)-trial_onsets(trialidx)+1;
end

nan_spiketimes_by_trial=find(isnan(spiketimes_by_trial(1,:)));
if(isempty(nan_spiketimes_by_trial))
    %do nothing
elseif(any(diff(nan_spiketimes_by_trial)>1))
    error('Found spiketimes not matched to any trial that did not occur before the first trial onset. Why?')
else
    spiketimes_by_trial(:,nan_spiketimes_by_trial)=[];
    spike_assignments(nan_spiketimes_by_trial)=[];
    spike_start_idx_by_run=spike_start_idx_by_run-nan_spiketimes_by_trial(end);
end
if (spike_start_idx_by_run~=1)
    error('spike_start_idx_by_run should be 1 but it''s not.')
end

if (~options.load_as_temp || options.force_compute_quality_measures)
    did_quality=true;
    %compute quality measures
    if(~isempty(best_channels_phy_merge))
        %give metrics a copy of spike_clusters where clusters that were
        %assigned the same channel and unit number are merged so that
        %metrics computed from the combined clusters are returned.
        spike_clusters = clusts;
        [~,i1,i2]=unique(double(best_channel)+double(unit_numbers_phy)/10);
        sorted_ids=ids_phy(i1);
        [~,si_sorted_ids]=sort(sorted_ids);
        for i=1:length(i2)
            spike_clusters(spike_clusters==ids_phy(i))=sorted_ids(i2(i));
        end
        [cgs, uQ, cR, isiV] = sqKilosort.computeAllMeasures([job.results_path,suffix],spike_clusters);
        mesure_cids=unique(spike_clusters);
    else
        [cgs, uQ, cR, isiV] = sqKilosort.computeAllMeasures([job.results_path,suffix],[],all_units);
        mesure_cids=ids;
    end
else
    did_quality=false;
end

if(0)
    usa=unique(spike_assignments);
    for i=1:length(usa)
        ns(i)=sum(spike_assignments==usa(i));
    end
    [~,si]=sort(uQ,1,'descend');
    fprintf(['\tID','\t  Group','\t   NSpikes',' \tiso_distance',' \tcont_rate*100',' \tisi_vio_rate*100\n'])
    disp([ids(si) cgs(si)' ns(si)' uQ(si) cR(si)*100,isiV(si)'*100])
end

spike_start_idx_by_run(end+1)=length(spike_assignments)+1;

Kilosort_load_completed_job_params=options;

if exist(merge_results_file,'file')
    Kilosort_load_completed_job_params.merge_file=merge_results_file;
else
    Kilosort_load_completed_job_params.merge_file='None';
end

% best_channel = channel for each unit
% best_channels = list of channels with any unit(s)

% remapping if KS2 removed electrodes from sort
if exist('rm_chans', 'var')
    for i=1:length(best_channel)
        n = sum(removed_inds <= best_channel(i));
        best_channel(i) = kept_inds(best_channel(i));
    end
end

best_channels=unique(best_channel);


% these are each a cell array, one entry per channel
% contains a matrix with info for each unit assigned to that channel
gSingleRawFields.isolation=cell(length(best_channels),1);
gSingleRawFields.phy_contamination_pct=cell(length(best_channels),1);
gSingleRawFields.unit_type=cell(length(best_channels),1);

dbopen

gSingleRawFields.unit_start = repmat({0},[length(best_channels),1]);
if(~options.load_as_temp)
    % purge existing database entries that metadata about these recordings
    rawids = [];
    for run_idx=1:length(job.runs)
        sql=['SELECT id,parmfile,bad',...
            ' FROM gDataRaw WHERE parmfile =''',job.runs{run_idx},''''];
        S=mysql(sql);
        rawids(end+1)=S.id;
    end
    if options.purge_rawids
        dblink_purge_rawids(rawids);     
    end
    % if append_units, figure out the maximum unit number that already exists
    % on each channel. new units added for these data will be assigned to
    % values starting at maxunit+1
    
    if options.append_units
        for ii=1:length(best_channels)
            sql=['SELECT max(unit) as maxunit FROM gSingleCell',...
                ' WHERE cellid like "',site,'%" AND channum=',num2str(best_channels(ii))];
            sdata = mysql(sql);
            if ~isempty(sdata.maxunit)
                gSingleRawFields.unit_start{ii}=sdata.maxunit;
            else
                gSingleRawFields.unit_start{ii}=0;
            end
        end   
    end

end


all_idx=true(size(spike_assignments));

for run_idx=1:length(job.runs)
    sql=['SELECT id,parmfile,bad',...
            ' FROM gDataRaw WHERE parmfile =''',job.runs{run_idx},''''];
    S=mysql(sql);
    if S.bad && ~options.load_as_temp
        fprintf('%s marked as bad, skipping\n',job.runs{run_idx});
        continue
    end
    all_idx(:)=true;
    all_idx(1:spike_start_idx_by_run(run_idx)-1)=false;
    all_idx(spike_start_idx_by_run(run_idx+1):end)=false;
    this_baphy_source=[job.runs_root,filesep,job.runs{run_idx}];
    spiketimes_this_trial=spiketimes_by_trial(:,all_idx);
    ind=0;
    s=cell(size(best_channels));
    ctypes=cell(size(best_channels));
    for i=1:length(best_channels)
        
        units=all_units(best_channel==best_channels(i));
        if(isempty(best_channels_phy_merge))
            ut=nan(size(units));
            ut(ismember(units,sus))=1;
            ut(ismember(units,mus))=2;
            un=(1:length(units))+options.unit_numbering_offset;
            [ut,si]=sort(ut);
            units=units(si);
        else
            ut=unit_types_phy(best_channel==best_channels(i));
            un=unit_numbers_phy(best_channel==best_channels(i));
            [~,si]=sortrows([ut,un]);
            ut=ut(si);
            un=un(si);
            units=units(si);
        end
        [uuns,unique_idx]=unique(un);
        ut=ut(unique_idx);
                
        %don't sort un, it's used for gettign the unit index only
        for j=1:length(uuns)
            ind=ind+1;
            
            % add offset to account for existing units that provide the
            % starting point for appended units
            final_ui=uuns(j) + gSingleRawFields.unit_start{i};
            
            un_idx=(un==uuns(j));
            idx=ismember(spike_assignments(all_idx),ids(units(un_idx)));
            %idx(spiketimes
            if ispc
                s{i}(final_ui,1).sorter=getenv('username');
            else
                s{i}(final_ui,1).sorter=getenv('USER');
            end
            s{i}(final_ui).primary=1;
            if ~options.delete_existing_spkfile && options.append_units
                s{i}(final_ui).primary=0;
                % If appending and not deleting, set primary to 0 to keep existing sort indexes
                % as they are, so database doen't need to be updated because old sorts got moved to sortidx2
            end
            s{i}(final_ui).comments='Sorted by KiloSort, manually curated by phy. Type is: ';
            if(~isempty(options.keep_Kids) && ~any(ismember(ids(units(un_idx)),options.keep_Kids)))
                ut(j)=2;
            end
            if(ut(j)==1)
                s{i}(final_ui).comments=[s{i}(final_ui).comments 'SU'];
                s{i}(final_ui).sortparameters.KiloSort_type='SU';
            else
                s{i}(final_ui).comments=[s{i}(final_ui).comments 'MU'];
                s{i}(final_ui).sortparameters.KiloSort_type='MU';
            end
            gSingleRawFields.unit_type{i}(final_ui)=ut(j);
            s{i}(final_ui).unitSpikes=spiketimes_this_trial(:,idx);
            s{i}(final_ui).Template=mean(templates(:,:,ismember(all_units,units(un_idx))),3);
            if isempty(mean_waveforms_phy)
                s{i}(final_ui).MeanWaveform=[];
            else
                s{i}(final_ui).MeanWaveform=mean_waveforms_phy(:,ismember(all_units,units(un_idx)));
            end
            s{i}(final_ui).env=[ones(size(s{i}(final_ui).Template)),repmat(s{i}(final_ui).Template,1,2)];
            s{i}(final_ui).Ncl=length(uuns) ;%clusters in this channel, redundant
            s{i}(final_ui).xaxis=[];
            s{i}(final_ui).sortparameters.SaveSorter=0;
            s{i}(final_ui).sortparameters.Kilosort_job_source=job_file;
            s{i}(final_ui).sortparameters.Kilosort_load_completed_job_params=Kilosort_load_completed_job_params;
            s{i}(final_ui).sortparameters.KiloSort_ids=ids(units(un_idx));
            s{i}(final_ui).mfilename=this_baphy_source;
            
            % save celltypes and associated stuff only for the appropriate 
            % channels... CRH
            if ~(isempty(fieldnames(wft_struct)))
                fns = fieldnames(wft_struct);
                for k = 1:length(fns)
                    key = char(fns(k));
                    ctypes{i}(final_ui).(key) = wft_struct.(key)(units(j));
                    if length(int2str(best_channels(i)))>1
                        ctypes{i}(final_ui).cellid= [site,'-',num2str(best_channels(i)),'-', num2str(final_ui)]; 
                    else
                        ctypes{i}(final_ui).cellid=[site,'-0',num2str(best_channels(i)),'-', num2str(final_ui)];
                    end

                end
            end
            
            if(did_quality)
                if(~isempty(best_channels_phy_merge))
                    metric_idx=find(ismember(mesure_cids,ids(units(un_idx))));
                    if(isempty(metric_idx))
                        error('Metrics were not computed for cluster id: %d. Why?',ids(units(un_idx)))
                    end
                    metric_idx(2:end)=[];
                else
                    metric_idx=units(un_idx);
                end
 
                % save the KS clusterID associated with this cellid
                % allClusterIDs = unique(clusts);
                allClusterIDs = mesure_cids;
                thisClusterID = allClusterIDs(metric_idx);
                s{i}(final_ui).ks_clusterID=thisClusterID;
                
                s{i}(final_ui).isolation_distance=uQ(metric_idx);
                s{i}(final_ui).contamination_perc=cR(metric_idx)*100;
                s{i}(final_ui).isi_violation_perc=isiV(metric_idx)*100;
                if(~isempty(snrs))
                    %iso perc is the error function of the maximum snr (across channels)
                    %for units that were assigned multiple clusters, be conservative and take the minumum
                    s{i}(final_ui).isolation_perc=min(100*erf(max(snrs(:,units(un_idx)),[],1)/2));
                else
                    s{i}(final_ui).isolation_perc=NaN;
                end
                
                % add ks cluster id to database for easier updating of db
                % tables later on
                gSingleRawFields.kilosort_cluster_id{i}(final_ui)=s{i}(final_ui).ks_clusterID;
                % likewise, add sorting path to reference the sorting job
                % for each cellid
                gSingleRawFields.ksJobPath{i}{final_ui}=job.results_path;
                
                
                gSingleRawFields.phy_contamination_pct{i}(final_ui)=s{i}(final_ui).contamination_perc;
                gSingleRawFields.isolation{i}(final_ui)=100-s{i}(final_ui).contamination_perc;
                if(ut(j)==1)
                    %if manually labeled as a SU, force isolation to be at least 95
                    gSingleRawFields.isolation{i}(final_ui)=max(gSingleRawFields.isolation{i}(final_ui),95);
                else
                    %if manually labeled as a MU, force isolation to be at most 85
                    gSingleRawFields.isolation{i}(final_ui)=min(gSingleRawFields.isolation{i}(final_ui),85);
                end
                s{i}(final_ui).isolation=gSingleRawFields.isolation{i}(final_ui);
                
                if(~isempty(options.keep_Kids) && ~any(ismember(ids(units(un_idx)),options.keep_Kids)))
                    gSingleRawFields.isolation{i}(final_ui)=5;
                    s{i}(final_ui).isolation=5;
                end
            else
                gSingleRawFields.phy_contamination_pct{i}(final_ui)=NaN;
                gSingleRawFields.isolation{i}(final_ui)=NaN;
                s{i}(final_ui).isolation_distance=NaN;
                s{i}(final_ui).contamination_perc=NaN;
                s{i}(final_ui).isi_violation_perc=NaN;
                s{i}(final_ui).isolation=NaN;
            end
            
            if(run_idx==1)
                if length(s{i}(final_ui).sortparameters.KiloSort_ids)==1
                    kstr=' %d';
                else
                    kstr=['s [',repmat('%d, ',1,length(s{i}(final_ui).sortparameters.KiloSort_ids))];
                    kstr=[kstr(1:end-2),']'];
                end
                fprintf(['Kilosort Id',kstr,' -> %s Channel %d, Unit %d \n'],s{i}(final_ui).sortparameters.KiloSort_ids,s{i}(final_ui).sortparameters.KiloSort_type,best_channels(i),final_ui)
            end
        end
    end
    
    if(options.load_as_temp)
        savfile=[fileparts(tempname) filesep job.runs{run_idx}(1:end-2) '.spk.mat'];
    else
        savfile = [job.runs_root filesep 'sorted' filesep job.runs{run_idx}(1:end-2) '.spk.mat'];
    end
    if(options.delete_existing_spkfile && exist(savfile,'file'))
        file_exists=true;
        while file_exists
            lastwarn('')
            delete(savfile)
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg) && strcmp(warnId,'MATLAB:DELETE:Permission')
                a=questdlg(['Warning! Permission denied when trying to delete ' savfile '. Delete manually usuing sudo then Continue. Or Cancel saving.'],'Permission Denied','Continue','Cancel','Cancel');
                if strcmp(a,'Cancel')
                    savfiles=[];
                    fprintf('\n Cancelled')
                    return
                end
            else
                file_exists = exist(savfile,'file');
            end
        end
    end
    extras.StimTagNames=exptparams_{run_idx}.TrialObject.ReferenceHandle.Names;
    extras.npoint=job.nSamplesBlocks(run_idx);
    extras.sweeps=exptparams_{run_idx}.TotalRepetitions;
    fprintf('savespikes.m: Saving spike data to: %s...\n',savfile)
    %format=6; %even if baphy format, use savespikes_do_save format 6 for saves from Kilosort/phy
    format=7; %format 7 means sort_idx is no longer used, length(sortinfo{chan}) should always be 1.
    
    switch job.runinfo(1).evpv
        case 5
            rate=job.fs;
        case 6
            rate=job.sampleRate(run_idx);
    end
    nrec=runs_per_trial(run_idx);
    extras.numChannels=job.Nchan;
    
    maxunits=1;
    for jj=1:length(s)
        maxunits=max(maxunits,size(s{jj},1));
    end
    
    goodtrials=repmat({''},max(cellfun(@length,s)),1); % OR '1:x' OR 'x:y' for excluding trials
    
    sortidxs=savespikes_do_save(savfile,s,best_channels,extras,format,rate,nrec);

    if(~options.load_as_temp)
        % need to talk to the celldb database
        fprintf('Saving to database\n')
        r=dbGetConParms();
        rtemp=r;
        rtemp.DB_SERVER='hyrax.ohsu.edu';
        
        % temporarily connect to NEMS db
        %dbSetConParms(rtemp);
        
        dblink_phy(savfile,s,best_channels,sortidxs,extras,gSingleRawFields,goodtrials);
        %matchcell2file_fn({savfile},{this_baphy_source},s,best_channels,extras,gSingleRawFields,goodtrials);
        
        % Now, save celltypes information to db
        % have to use cellids here (maybe there is a better way to do
        % this? integrate with matchcell2file...?
        if ~(isempty(fieldnames(wft_struct)))
            for kk = 1:length(best_channels)
                if length(ctypes{kk})>1
                    for c = 1:length(ctypes{kk})
                        cellid = ctypes{kk}(c).cellid;
                        t = rmfield(ctypes{kk}(c), 'cellid');
                        if length(cellid) ~= 0
                            dbWriteMetaData(cellid,t)
                        end
                    end
                else
                    cellid = ctypes{kk}(1).cellid;
                    t = rmfield(ctypes{kk}(1), 'cellid');
                    dbWriteMetaData(cellid,t)
                end
            end       
        end
        
        % restore old db connection
        %dbSetConParms(r);
    end
    
    %copy over additional fields
    seplocs = strfind(this_baphy_source,filesep);
    fname = [this_baphy_source(seplocs(end-1)+1:end) '.evp'];
    exptevents=exptevents_{run_idx};
    save(savfile,'-append','exptevents','fname')
    savfiles{run_idx}=savfile;
end

if(~options.load_as_temp && isempty(strfind(job_file,'completed')))
    dest=strrep(job_file,'in_progress','completed');
    if ~(exist(dest)==2)
        UTmkdir(dest)
        movefile(job_file,dest)
    end
end

end
