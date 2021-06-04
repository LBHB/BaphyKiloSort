rez=load([job.results_path_temp,'/rez.mat'],'dshift','ops','F','F0');
batch_starts=1:rez.ops.NT:rez.ops.tend;
batch_mids=batch_starts+rez.ops.NT/2;
run_starts = [1 cumsum(rez.ops.nSamplesBlocks(1:end-1))];
figure;plot(batch_starts/rez.ops.fs/60,fastsmooth(rez.dshift,3,2,1))
hold on;
yl=get(gca,'YLim');
for i=1:length(run_starts)
    plot(run_starts([i i])/rez.ops.fs/60,yl,'--','Color',[.5 .5 .5]);
    str=[rez.ops.runs{i}(8:9),'\_',rez.ops.runs{i}(13:15)];
    th(i)=text(run_starts(i)/rez.ops.fs/60,yl(2),str,'HorizontalAlignment','left','VerticalAlignment','top');
end
utitle(job.name)
ylabel('drift (um)')
xlabel('Time (min)')
set(gca,'Box','off','TickDir','out')

uiopen([job.results_path filesep 'dift_map.fig'],1)
set(findall(gcf,'LineStyle','--'),'Visible','off')

if 0 
    i=1;
nf=max(rez.F(:))/255;
imN=rez.F/nf;
imN0=rez.F0/nf;
figure;imagesc(imN0);

figure;im=image(imN(:,:,i));
title(num2str(i))

while(1)
    k = waitforbuttonpress;
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    value = double(get(gcf,'CurrentCharacter'));
    if value == 28
        i=i-1;
    elseif value == 29
        i=i+1;
    else
        a=inputdlg('Frame #');
        i=str2double(a{1});
    end
    set(im,'CData',imN(:,:,i))
    title(num2str(i))
end
    
end