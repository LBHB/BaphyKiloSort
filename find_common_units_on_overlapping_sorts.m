function find_common_units_on_overlapping_sorts(job_pth1,job_pth2)
%find_common_units_on_overlapping_sorts(job_pth1,job_pth2)
%
% Finds matching clusters between two jobs that both include the same run
% Used to find matches in a short job that is a subset of a longer job
% 
% For each cluster in job1, computes the trial PSTH for that cluster and
%  finds the cluster in job2 that has the highest cross-correlation (only among labeled job2 clusters)
%  Creates 3 TSV files that will be used by phy to show the matches
%    cluster_XC - The cross-correlation value of the match
%    cluster_matchID - The KiloSort id for the match (in job2)
%    cluster_matchLabel - The label of the match (in job2)
%
% Inputs:
%   job_pth1 - The full path to the short job
%   job_pth2 - The full path to the long (reference) job
%   You can easily get these by pressing 'load job to workspace' in KiloSort browser

 
fs=20; %Sampling rate for PSTHs over which x-corrs are calculated

job1=load(job_pth1);
job2=load(job_pth2);
%rez = load([job1.results_path_temp filesep 'rez.mat'])
spike_times1=double(readNPY([job1.results_path filesep 'spike_times.npy']));
spike_times2=double(readNPY([job2.results_path filesep 'spike_times.npy']));
spike_clusters1=readNPY([job1.results_path filesep 'spike_clusters.npy']);
spike_clusters2=readNPY([job2.results_path filesep 'spike_clusters.npy']);
ids1 = unique(spike_clusters1);

cluster_group_file2=[job2.results_path filesep 'cluster_group.tsv'];
fid= fopen(cluster_group_file2);
C = textscan(fid,'%d %s',-1,'Headerlines',1,'Delimiter','\t');
fclose(fid);
ids2=C{1};
group2=C{2};
% ids2 = unique(spike_clusters2);
% group2=repmat({''},size(ids2));
% for i=1:length(ids2_)
%     group2{ids2_(i)==ids2}=group2_{i};
% end
% Keep only group == 'good' | group == 'mua' and update ids
% su_mu_mask = cellfun(@(x) (strcmp(x,'good')|strcmp(x, 'mua')), group2);
% group2 = group2(su_mu_mask);
% ids2 = ids2(su_mu_mask);


run_start_idx = find(strcmp(job1.runs{1},job2.runs));
spike_time_offset = job2.trial_onsets_{run_start_idx}(1) - job1.trial_onsets_{1}(1);
max_time = max(spike_times1);
binsize = job1.sampleRate/fs;
if length(binsize)>1
    binsize = unique(binsize);
    warning("SP(7/8/25): multiple binsizes, choosing unique binsize: %d\n", binsize)
end
tbin_edges = 0:binsize:(max_time+binsize); % vector of time bin edges (for histogram)

job1_raster = zeros(length(tbin_edges),length(ids1));
for i=1:length(ids1)
    job1_raster(:,i) = histc(spike_times1(spike_clusters1==ids1(i)),tbin_edges); % get spike counts for each bin
end

job2_raster = zeros(length(tbin_edges),length(ids2));
for i=1:length(ids2)
    job2_raster(:,i) = histc(spike_times2(spike_clusters2==ids2(i)) - double(spike_time_offset),tbin_edges); % get spike counts for each bin
end

xc=zeros(length(ids1),length(ids2));
for i=1:length(ids1)
    for j=1:length(ids2)
       xct=corrcoef(job1_raster(:,i),job2_raster(:,j)); 
       xc(i,j)=xct(1,2);
    end
end

for i=1:length(ids1)
    [best_corr(i),mi]=max(xc(i,:));
    best_match(i) = ids2(mi);
    best_matchLabel{i} = group2{mi};
end

fi_corr = fopen([job1.results_path filesep 'cluster_XC.tsv'],'w');
fprintf(fi_corr, 'cluster_id%sXC', char(9));
fprintf(fi_corr, char([13 10]));

fi_match = fopen([job1.results_path filesep 'cluster_matchID.tsv'],'w');
fprintf(fi_match, 'cluster_id%smatchID', char(9));
fprintf(fi_match, char([13 10]));

fi_matchLabel = fopen([job1.results_path filesep 'cluster_matchLabel.tsv'],'w');
fprintf(fi_matchLabel, 'cluster_id%smatchLabel', char(9));
fprintf(fi_matchLabel, char([13 10]));

for i=1:length(ids1)
    fprintf(fi_corr, '%d%s%.2f', ids1(i), char(9),best_corr(i));
    fprintf(fi_corr, char([13 10]));
    
    fprintf(fi_match, '%d%s%d', ids1(i), char(9),best_match(i));
    fprintf(fi_match, char([13 10]));
    
    fprintf(fi_matchLabel, '%d%s%s', ids1(i), char(9),best_matchLabel{i});
    fprintf(fi_matchLabel, char([13 10]));    
end
fclose(fi_corr);
fclose(fi_match);
fclose(fi_matchLabel);
[had_error,msg]=unix(['chmod a+rw ' job1.results_path filesep 'cluster_XC.tsv']);
if had_error
    error(msg)
end
[had_error,msg]=unix(['chmod a+rw ' job1.results_path filesep 'cluster_matchID.tsv']);
if had_error
    error(msg)
end
[had_error,msg]=unix(['chmod a+rw ' job1.results_path filesep 'cluster_matchLabel.tsv']);
if had_error
    error(msg)
end