function ops = concatenateOpenEphysBinary(ops)

force_rerun = false;

% generates a cell (channel) of structs (blocks) of .continuous files across
% all different directories corresponding to different experimental blocks
ch=load(ops.chanMap);
chans=sort(ch.chanMap);
fs=cell(length(ops.root),1);

for j = 1:ops.NchanTOT
    ops.chanMap_KiloRaw(ch.chanMap==chans(j))=j;
end

for k = 1:length(ops.root)
    % finds the AP recording node file
    bin_path_parts = strsplit(ops.runinfo(k).json_file, filesep);
    fs{k} = fullfile(ops.root{k}, bin_path_parts{end-3:end-1}, ...
                       'continuous', 'Neuropix-PXI-100.0', 'continuous.dat');
                   
    %calculates the number of samples for each file
    s = dir(fs{k});  
    nsamps = s.bytes / 2 / ops.NchanTOT;
    ops.nSamplesBlocks(1, k) = nsamps;
end

nBlocks     = length(ops.root);
nSamples    = 1024;  % fixed to 1024 for now!

% input .dat files is in data base temp
% output .bin file in the same folder.
fname = ops.fbinary;

% check if a reasonable binary file already exists
if exist(fname, 'file') == 2 && ~force_rerun
    %ensures the number of samples of the concatenated binary matches the
    %source binaries
    in_nsamps = sum([ops.nSamplesBlocks(1,:)]);
    out_nsamps = dir(fname).bytes / 2 / ops.NchanTOT;
    assert(length(in_nsamps) == length(out_nsamps));    
    
    fprintf('concatenated file found, loading.\n')
    
else
    UTmkdir(fname);

    fidout      = fopen(fname, 'w');
    if(fidout==-1)
        error(['Could not open file: ',fname])
    end



    fprintf('Concatenating Open Ephys binary files to a single local binary file.\n')
    tic
    for k = 1:nBlocks
        fprintf('file %d of %d \n', k, nBlocks')

        fid = fopen(fs{k});

        while ~feof(fid)

            % N channels x M samples 16-bit integers in little-endian format. 
            % Data is saved as ch1_samp1, ch2_samp1, ... chN_samp1, 
            % ch1_samp2, ch2_samp2, ..., chN_sampM. 
            samples = fread(fid, [ops.NchanTOT, nSamples * 1000],'int16', 'ieee-le');

             if(isfield(ops,'common_rejection_mode'))
                     switch ops.common_rejection_mode
                         case 'none'
                         case 'mean'
                             samples=samples-repmat(mean(samples),size(samples,1),1);
                         case 'median'
                             samples=samples-repmat(median(samples),size(samples,1),1);
                         otherwise
                             error(['Unknown comon rejection mode ',ops.common_rejection_mode])
                     end
             end
            fwrite(fidout, samples, 'int16');

        end

        fclose(fid); 
    end

    fclose(fidout);
    fprintf('Done concatenating, ')
    toc
    fprintf('\n')
    
end
end
