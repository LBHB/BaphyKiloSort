function parent_handle = Kilosort_browser(parent_handle)

global SERVER_PATH 
KSpaths=UTkilosort_paths();
h = 610;  % Total window height
w = 950; % Total window width
ht=500;
pad=3;
ts=20;
do_buttons_hpos=300;
do_buttons_hpos2=550;
show_completed=1;
show_incomplete=1;
if ~exist('parent_handle', 'var')
    parent_handle = figure('Menubar','none', 'Resize','off', ...
        'Units','pixels', 'Position', [20 50 w h],...
        'Name', 'Kilosort Jobs Browser', 'NumberTitle', 'off',...
        'HandleVisibility', 'off');
end

if isempty(which('TableSorter'))
    TPATH = fileparts(which('narf_set_path'));
    javaaddpath([TPATH '/libs/TableSorter.jar']);
end
if isempty(which('TableColumnAdjuster'))
    TPATH = fileparts(which('narf_set_path'));
    javaaddpath([TPATH '/libs/TableColumnAdjuster.jar']);
end

uicontrol('Parent', parent_handle, 'Style', 'text', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'Name:', ...
    'BackgroundColor',[.8 .8 .8],...
    'Position', [pad ht+4*ts 50 ts]);
handles.name = uicontrol('Parent', parent_handle, 'Style', 'edit', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', '', ...
    'Max', 1, 'Min', 1, ...
    'Position', [2*pad+50 ht+4*ts 200 ts], ...
    'Callback', @reload_job_table);

handles.view_in_phy_button = uicontrol('Parent', parent_handle,...
    'Style', 'pushbutton', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'View in Phy', ...
    'Max', 1, 'Min', 1, ...
    'Position', [do_buttons_hpos ht+4*ts 200 ts], ...
    'Callback', @view_in_phy_callback,'BUTTONDOWNFCN',@view_in_phy_callback);
handles.save_to_db_button = uicontrol('Parent', parent_handle,...
    'Style', 'pushbutton', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'Save to database', ...
    'Max', 1, 'Min', 1, ...
    'Position', [do_buttons_hpos ht+3*ts 200 ts], ...
    'Callback', @save_to_db_callback);

handles.delete_job_button = uicontrol('Parent', parent_handle,...
    'Style', 'pushbutton', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'Delete', ...
    'Max', 1, 'Min', 1, ...
    'Position', [do_buttons_hpos2 ht+4*ts 200 ts], ...
    'Callback', @delete_job_callback);

handles.show_log_button = uicontrol('Parent', parent_handle,...
    'Style', 'pushbutton', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'Show Log', ...
    'Max', 1, 'Min', 1, ...
    'Position', [do_buttons_hpos2 ht+3*ts 200 ts], ...
    'Callback', @show_log_job_callback);

handles.load_job_button = uicontrol('Parent', parent_handle,...
    'Style', 'pushbutton', 'Units', 'pixels',...
    'HorizontalAlignment', 'left', 'String', 'Load job to workspace', ...
    'Max', 1, 'Min', 1, ...
    'Position', [do_buttons_hpos2 ht+2*ts 200 ts], ...
    'Callback', @load_job_callback);

handles.job_table = uitable('Parent', parent_handle, ...
    'Enable', 'on', 'Units', 'pixels',...
    'RowName', [], ...
    'ColumnWidth', {700,80,80,80}, ...
    'ColumnName', {'Name', 'Sorted', 'Saved in Phy','Saved to|database'}, ...
    'Position', [pad pad w-pad*2 ht]);


% Configure the job table selection to update the center and right panels
hJS = findjobj(handles.job_table);
hJT = hJS.getViewport.getView;
hJT.setNonContiguousCellSelection(false);
hJT.setColumnSelectionAllowed(false);
hJT.setRowSelectionAllowed(true);
hJT.setSelectionMode(0); % Allow only a single row to be selected at once
hJTcb = handle(hJT, 'CallbackProperties');
set(hJTcb, 'MouseReleasedCallback', {@job_table_row_selected, parent_handle});
set(hJTcb, 'KeyPressedCallback', {@job_table_row_selected, parent_handle});
jobs=[];
ji=[];
%reload_job_table

    function job_table_row_selected(Hobj,ev,a)
        ji = Hobj.getSelectedRows()+1;
        j=load([jobs(ji).pth filesep jobs(ji).name],'results_path');
        r_pth=j.results_path;
        cmd=['phy_lbhb template-gui ' r_pth filesep 'params.py\n'];
        fprintf(cmd)
    end

    function reload_job_table(Hobj,ev)
        enable_buttons(0)
        jobs=[];
        name_filter=get(handles.name,'String');
        if(show_incomplete)
            pth=[KSpaths.jobs_path filesep 'in_progress'];
            jobs=dir(pth);jobs(1:2)=[];
            jobs(arrayfun(@(x)length(x.name)<5,jobs))=[];
            jobs(arrayfun(@(x)~strcmp(x.name(end-3:end),'.mat'),jobs))=[];
            if(~isempty(name_filter))
                jobs(arrayfun(@(x)isempty(strfind(x.name,name_filter)),jobs))=[];
            end
            if(~isempty(jobs))
                for i=1:length(jobs)
                    tmp=load([pth filesep jobs(i).name],'status','results_path');
                    jobs(i).sorted=logical(tmp.status);
                    jobs(i).saved_in_phy=exist([tmp.results_path filesep 'cluster_group.tsv'],'file')==2;
                    ss=strsep(strrep(jobs(i).name,'.mat',''),'_');
                    jobs(i).pen=ss{1};
                    jobs(i).rn1=ss{2};
                    if ~isnumeric(jobs(i).rn1)
                       ind=length(ss{2})-1;
                       while ind>0
                           jobs(i).rn1=str2double(ss{2}(1:ind));
                           if ~isnan(jobs(i).rn1)
                               break
                           end
                           ind=ind-1;
                       end
                       if ind==0
                           jobs(i).rn1=NaN;
                           warning(['Couldn''t figure out run number for ',jobs(i).name,'. Setting to NaN.'])
                       end
                    end
                end
                sf=cell(length(jobs),2);
                [sf{:,1}]=deal(jobs.pen);
                [sf{:,2}]=deal(jobs.rn1);
                [~,si]=sortrows(sf,[1 2]);
                jobs=jobs(si);
                [jobs.saved_in_db]=deal(false);
                [jobs.pth]=deal(pth);
            else
               jobs=[]
            end
        end
        if(show_completed)
            pth=[KSpaths.jobs_path filesep 'completed'];
            Cjobs=dir(pth);Cjobs(1:2)=[];
            Cjobs(arrayfun(@(x)length(x.name)<5,Cjobs))=[];
            Cjobs(arrayfun(@(x)~strcmp(x.name(end-3:end),'.mat'),Cjobs))=[];
            if(~isempty(name_filter))
                Cjobs(arrayfun(@(x)isempty(strfind(x.name,name_filter)),Cjobs))=[];
            end
            if(~isempty(Cjobs))
                for i=1:length(Cjobs)
                    tmp=load([pth filesep Cjobs(i).name],'status','results_path');
                    Cjobs(i).sorted=logical(tmp.status);
                    Cjobs(i).saved_in_phy=exist([tmp.results_path filesep 'cluster_group.tsv'],'file')==2;
                    ss=strsep(Cjobs(i).name,'_');
                    Cjobs(i).pen=ss{1};
                    Cjobs(i).rn1=ss{2};
                end
                sf=cell(length(Cjobs),2);
                [sf{:,1}]=deal(Cjobs.pen);
                [sf{:,2}]=deal(Cjobs.rn1);
                [~,si]=sortrows(sf,[1 2]);
                Cjobs=Cjobs(si);
                [Cjobs.saved_in_db]=deal(true);
                [Cjobs.pth]=deal(pth);
                jobs=[jobs;Cjobs];
            end
        end
        
        
        l = length(jobs);
        c = cell(l,4);
        for i = 1:l
            c{i,1} = jobs(i).name(1:end-4);
            c{i,2} = jobs(i).sorted;
            c{i,3} = jobs(i).saved_in_phy;
            c{i,4} = jobs(i).saved_in_db;
        end
        set(handles.job_table, 'Data', c);
        a=2;
        enable_buttons(1)
    end

    function view_in_phy_callback(Hobj,ev)
        if(isempty(ji))
            return
        end
        st=get(get(Hobj,'Parent'),'SelectionType');
        enable_buttons(0);
        j=load([jobs(ji).pth filesep jobs(ji).name]);
        r_pth=j.results_path;
        cmd=['phy_lbhb template-gui ' r_pth filesep 'params.py'];
        if(strcmp(st,'alt'))
            cmd=[cmd,' &'];
        end
        fprintf(['Loading using command: ',cmd,'\n'])
        err_type=0;
        try
            [s,w]=unix(cmd,'-echo');
        catch err
            err_type=1;
        end
        if err_type==0 && ~strcmp(st,'alt')
            try
            %if phy was run blocking matlab, change permissions of new
            %files phy created so other users can look at this file too.
            [w2,s]=unix(['chmod -R 777 ',r_pth]); if w2, error(s), end
            catch err
                err_type=2;
            end
        end
        
        enable_buttons(1);
        if err_type==1
            rethrow(err);
        elseif err_type==2
            warning('Could not set permissions to be open for all files in %s. This is probably because someone else owns some files here.',r_pth)
        end
    end
    function save_to_db_callback(Hobj,ev)
        %enable_buttons(0);
        [s, button] = settingsdlg(...
            'Description', ['Saving ',jobs(ji).name,...
            '. Set save parameters:'],...
            'title'      , 'Kilosort save options',...
            'WindowWidth',650,...
            'ControlWidth',350,...
            'separator','Saving to Database',...
            {'Purge rawids (clear celldB database)'; 'purge_rawids'}, true,...
            {'Delete existing spk.mat file(s)'; 'delete_existing_spkfile'}, true,...
            {'Append unit numbers (for each channel, create a new unit number if a unit already exists for any run)'; 'append_units'}, true,...
            {'Don''t save these ids (don''t save these units)';'remove_ids'},'',...
            'separator','Temporary Saving',...
            {'Save to temp (save spk.mat file in /tmp/, don''t write database)'; 'load_as_temp'}, false,...
            {'Force compute quality (Computre quality metrics, by default not computed if saving to temp)'; 'force_compute_quality_measures'}, false);
        
        if(strcmp(button,'OK'))
            savfiles=UTkilosort_load_completed_job([jobs(ji).pth filesep jobs(ji).name],...
                0, s.load_as_temp, s.force_compute_quality_measures, ...
                s.append_units,s.delete_existing_spkfile, s.purge_rawids, s.remove_ids);
            reload_job_table();
        end
        
    end
    function enable_buttons(en)
        list={'name','view_in_phy_button','save_to_db_button',...
            'job_table','delete_job_button'};
        for i=1:length(list)
            if(en)
                set(handles.(list{i}),'Enable','on')
            else
                set(handles.(list{i}),'Enable','off')
            end
        end
        pause(.1)
    end
    function show_log_job_callback(Hobj,ev)
        if(isempty(ji))
            return
        end
        j=load([jobs(ji).pth filesep jobs(ji).name],'queueidx');
        if isempty(fieldnames(j))
            warning(['queueidx not found. These were only more recently ',...
                'stored in the job file so they''re not available for ',...
                'older files.']);
        else
            fn=sprintf(['%sweb',filesep,'celldb',filesep,'queue',filesep,...
            '%d',filesep,'%d.out'],SERVER_PATH,floor(j.queueidx/1000)*1000,...
            j.queueidx);
            if exist(fn,'file')
                fprintf(sprintf('-------------- Log for %s -----------------\n',...
                [jobs(ji).pth filesep jobs(ji).name]))
                type(fn)
                fprintf('-------------- End Log -----------------\n')
            else
                warning(sprintf('Log file %s does not exist',fn))
            end
        end
        
    end
    function load_job_callback(Hobj,ev)
        if(isempty(ji))
            return
        end
        j=load([jobs(ji).pth filesep jobs(ji).name]);
        assignin('base','job',j)
        fprintf(sprintf('Loaded %s to ''job''\n',...
                [jobs(ji).pth filesep jobs(ji).name]))
    end
    function delete_job_callback(Hobj,ev)
        if(isempty(ji))
            return
        end
        [s, button] = settingsdlg(...
            'Description', ['Deleting ',jobs(ji).name],...
            'title'      , 'Delete job',...
            {'Delete phy file'; 'delete_phy'},true,...
            {'Delete job file'; 'delete_mat'},true);           
        
        if(strcmp(button,'OK'))
            if(s.delete_phy)
                tmp=load([jobs(ji).pth filesep jobs(ji).name],'results_path');
                sl=strfind(tmp.results_path,filesep);
                [status,m,mid]=rmdir(tmp.results_path(1:sl(end)),'s');
                if(~status)
                    switch mid
                        case 'MATLAB:RMDIR:NotADirectory'
                            warning(m)
                        otherwise
                            error(m)
                    end
                end
            end
            if(s.delete_mat)
                delete([jobs(ji).pth filesep jobs(ji).name]);
            end
            reload_job_table();
        end 
        
end
end