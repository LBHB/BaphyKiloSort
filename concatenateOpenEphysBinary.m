function ops = concatenateOpenEphysBinary(ops)

force_rerun = false;

% generates a cell (channel) of structs (blocks) of .continuous files across
% all different directories corresponding to different experimental blocks
if isstruct(ops.chanMap)
    ch=ops.chanMap;
    ops.chanMap_KiloRaw = ops.chanMap.chanMap * 0;
    connected = logical(ops.chanMap.connected);

    ops.chanMap_KiloRaw(connected) = 1:ops.NchanTOT;
else
    ch=load(ops.chanMap);
    chans=sort(ch.chanMap);
    fs=cell(length(ops.root),1);
    
    for j = 1:ops.NchanTOT
        ops.chanMap_KiloRaw(ch.chanMap==chans(j))=j;
    end
end

ops.first_timestamps=zeros(1,length(ops.root));
for file_num = 1:length(ops.root)
    % finds the AP recording node file
    bin_path_parts = strsplit(ops.runinfo(file_num).json_file, filesep);
    % data_name depends on whether this is Neuropixels or Rhythm FPGA
    data_name=ops.runinfo(file_num).data_name;
    tsf = fullfile(ops.root{file_num}, bin_path_parts{end-3:end-1}, ...
                     'continuous', data_name, 'sample_numbers.npy');
    sample_numbers = readNPY(tsf);
    ops.first_timestamps(file_num) = sample_numbers(1);
    fs{file_num} = fullfile(ops.root{file_num}, bin_path_parts{end-3:end-1}, ...
                     'continuous', data_name, 'continuous.dat');
    fprintf('%d: %s t0=%d\n', file_num, fs{file_num}, ops.first_timestamps(file_num));
    %calculates the number of samples for each file
    s = dir(fs{file_num});
    nsamps = s.bytes / 2 / ops.NchanTOT;
    ops.nSamplesBlocks(1, file_num) = nsamps;
end

nBlocks     = length(ops.root);
% nSamples    = 1024;  % fixed to 1024 for now!
nSamples    = 256;  % SP changed so that filtfilt inside remove_mic_noise can be on gpu 

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
    for file_num = 1:nBlocks
        fprintf('file %d of %d \n', file_num, nBlocks')

        fid = fopen(fs{file_num});
        seg_num=0;
        while ~feof(fid)
            
            %rm_noise= 1;
            % N channels x M samples 16-bit integers in little-endian format. 
            % Data is saved as ch1_samp1, ch2_samp1, ... chN_samp1, 
            % ch1_samp2, ch2_samp2, ..., chN_sampM. 
            samples = fread(fid, [ops.NchanTOT, nSamples * 1000],'int16', 'ieee-le');
            % remove_mic_noise has been run for REI050, REI055, REI058
            if ops.rm_noise
                fprintf('Selected option to remove noise.\n')
                % if contains(ops.name, 'REI059')
                addpath('/auto/users/lbhb/Code/SP_matlab/');
                % fs=30e3;
                do_plot= 1;
                samples= remove_mic_noise(samples, ops.results_path_temp, ops.fs, do_plot, file_num, seg_num);
            end

            
            seg_num=seg_num+1;
            s = std(samples(:));
            %[seg_num, s]
            
            if s>350
                fprintf('cc=%d std %.3f over threshold, zeroing\n', seg_num, s);
                samples(:)=0;
                %keyboard
            end

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
