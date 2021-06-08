function UTkilosort2_plot_drift(rez,do_save)
% UTkilosort2_plot_drift(job) where job is the path to the job file
% UTkilosort2_plot_drift(rez) where rez is the loaded rez structure
% UTkilosort2_plot_drift(job,do_save) do_save is 0 by default
% UTkilosort2_plot_drift(rez,do_save)

if nargin==1
    do_save=0;
end
if ischar(rez)
    job=load(rez);
    rez=load([job.results_path_temp filesep 'rez.mat'],'ccb','ccbsort','ops');
end

figure('Position',[300 478 1000 420]);
subplot(1,2,1)
imagesc(rez.ccb, [-5 5]); drawnow
axis square
xlabel('batches')
ylabel('batches')
title('batch to batch distance')

subplot(1,2,2)
imagesc(rez.ccbsort, [-5 5]); drawnow
axis square
xlabel('sorted batches')
ylabel('sorted batches')
title('AFTER sorting')

ax=gca;
p=get(ax,'Position');
cb=colorbar(ax);
set(ax,'Position',p)
ylabel(cb,'Similarity')
set(gcf,'Color','w')

try
    set(gcf,'Name',job.name)
end

if do_save
    export_fig([rez.ops.results_path filesep 'drift.png'],'-a1','-m3','-q100');
end