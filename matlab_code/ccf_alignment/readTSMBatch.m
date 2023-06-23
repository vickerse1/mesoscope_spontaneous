function readTSMBatch(FileIN, FileOUT, XSize, YSize,TotalFrames,num_bins)
%% Adapted from readTSMFile by E. Abe on 06 Feb 2019




frames_per = TotalFrames/num_bins; 

%fseek(fid, 2880+(819200*0), 'bof');

for i  = 1:num_bins
    tic
    i
    fid = fopen(FileIN);
    fseek(fid, 2880+(819200*((i-1)*frames_per)), 'bof');
    this_batch = (char(num2str(i)));
    videoWF = fread(fid, [XSize*YSize, frames_per], '*int16');
    
    videoWF = reshape(videoWF,XSize,YSize,frames_per);
    videoWF =  (videoWF(:,2:end-1,:)); % Take out first and last column as they are extra data.
    savepath = [FileOUT '_' this_batch '.mat'];
    save(savepath,'videoWF','-v7.3');
    clear videoWF
    fclose(fid);
    
    toc
end

end 