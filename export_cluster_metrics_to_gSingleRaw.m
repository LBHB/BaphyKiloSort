% Simple script to load clusters / KS job and export spike isolation
% metrics for this job to celldb. 

% For, now this just means populating with isi_violation_pct. Might want to
% add other metrics at some point in the future. isi_violation now gets
% saved by default with new sorts

global SERVER_PATH

% CRH 08.16.2020
animal = 'Cordyceps';
site = 'CRD002a';
penetration = site(1:end-1);
run_nums = [15 16 17 18 19];
job_folder = 'CRD002a_15_16_17_18_19_KiloSort2_minfr_goodchannels_to0';
automerge = false;

if automerge
    results_path = [SERVER_PATH, 'daq', filesep, animal, filesep, ...
                        penetration, filesep, 'tmp', filesep, 'KiloSort',...
                        filesep, job_folder, filesep, ...
                        'results_after_automerge'];  
else
    results_path = [SERVER_PATH, 'daq', filesep, animal, filesep, ...
                        penetration, filesep, 'tmp', filesep, 'KiloSort',...
                        filesep, job_folder, filesep, ...
                        'results'];
end

merge_results_file=[results_path,filesep,'cluster_names.ts'];
if exist(merge_results_file,'file')
    clusts=readNPY([results_path, filesep, 'spike_clusters.npy']);
    error('Not implemented for merge cases yet')
else
    [cgs, uQ, cR, isiV, clusterID] = sqKilosort.computeAllMeasures(results_path);
    clusterID = clusterID - 1;
end

% get cellids corresponding to this kilosort run (i.e. corresponding to
% these run numbers
dbopen;
sql = ['SELECT DISTINCT gSingleRaw.cellid, sCellFile.respfile, gSingleRaw.kilosort_cluster_id, gSingleRaw.rawid FROM gSingleRaw', ...
            ' JOIN sCellFile ON (sCellFile.rawid=gSingleRaw.rawid) ', ...
            sprintf('WHERE sCellFile.cellid like %s', ['"%', site, '%"'])];
d = mysql(sql);
keep_rows = zeros(length(d), 1);
for ii = 1:length(d)
    sp1 = strsplit(d(ii).respfile, site);
    sp2 = strsplit(sp1{2}, '_');
    if ismember(str2double(sp2{1}), run_nums)
        keep_rows(ii) = 1;
    else
        keep_rows(ii) = 0;
    end
end
d = d(keep_rows==1);
       

% Now, update table by adding ISI violations for each cluster
for ii = 1:length(d)
   kid = d(ii).kilosort_cluster_id;
   if ~isempty(kid)
       isi_vio = isiV(clusterID==kid) * 100;
       sql = [sprintf('UPDATE gSingleRaw SET isi_contamination_pct=%f WHERE ', isi_vio), ...
                    sprintf('cellid="%s" and kilosort_cluster_id=%d', d(ii).cellid, kid)];
       fprintf('Updating isi_contamination_pct for cellid: %s, rawid: %d', d(ii).cellid, d(ii).rawid)
       fprintf('\n')
       mysql(sql)
   end
    
end




