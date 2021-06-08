function [template, templateMinIndex,was_split,was_merged] = UTkilosort_get_template(clu,template_ids,clusterID,rez)
% returns the spatiotemporal template that is the mean of
% all templates used for a given cluster.  Also returns which channel
% of that template has the maximally minimum value.


templatesUsed = unique(template_ids(clu==clusterID));
%   load(fullfile(penDir,'rez.mat'),'Wraw'); %the raw templates are in rez.Wraw
was_merged=length(templatesUsed)>1;
if(~isequal(clu==clusterID,ismember(template_ids,templatesUsed)))
    warning(['Cluster ID: ',num2str(clusterID),' was split. Template uncertain.'])
    was_split=true;
else
    was_split=false;
end
template=mean(rez.Wraw(:,:,templatesUsed+1),3);

%Which channel had the maximally minimum spike
[~,templateMinIndex] = min(min(template,[],2));
end