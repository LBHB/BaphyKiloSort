% function queueidx=dbaddqueue(rundataid,progname,parmstring,[allowqueuemaster],[note],[user]);
%
% add an entry to tQueue
%
% Created SVD 6/01
%
function queueidx=dbaddqueue(rundataid,progname,parmstring,allowqueuemaster,note,user,GPU_job);

dbopen;

global DBUSER
if ~exist('allowqueuemaster','var'),
   allowqueuemaster=0;
end
if ~exist('note','var'),
   note='';
end
if ~exist('user','var'),
   user=DBUSER;
end

if ~exist('GPU_job','var'),
    GPU_job=0;
end
y=year(now);
mon=month(now);
d=day(now);
h=hour(now);
min=minute(now);
s=round(second(now));
queuedate=sprintf('%d-%.2d-%.2d %.2d:%.2d:%.2d',y,mon,d,h,min,s);

[aff,queueidx]=sqlinsert('tQueue',...
                         'rundataid',rundataid,...
                         'progname',progname,...
                         'parmstring',parmstring,...
                         'priority',1,...
                         'queuedate',queuedate,...
                         'allowqueuemaster',allowqueuemaster,...
                         'note',note,...
                         'GPU_job',GPU_job,...
                         'user',user);
   
%fprintf('Created queue entry %d.\n',queueidx);
dbevent(6,queueidx);
