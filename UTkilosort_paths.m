function ops=UTkilosort_paths
global SERVER_PATH
if(isempty(SERVER_PATH))
    %LAS: sometimes SERVER_PATH gets emptied. Didn't figure out why. This fixes it.
    baphy_set_path
    if(isempty(SERVER_PATH))
        error('SERVER_PATH is empty!')
    end
end
ops.jobs_path=[SERVER_PATH 'code' filesep 'KiloSort' filesep 'jobs'];