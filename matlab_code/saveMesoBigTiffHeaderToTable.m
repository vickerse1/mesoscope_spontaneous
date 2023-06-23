prompt='Enter path: '
path=input(prompt,'s');
path=strcat(path,'\')

prompt='Enter filename: '
filename=input(prompt,'s');

file=strcat(path,filename);

[header,Aout,imgInfo,rawStream] = scanimage.util.opentif(file);

%frameTimestamps_sec=header.frameTimestamps_sec;

starttime=string(datetime(header.epoch{1}))

%target_file_csv=strcat(file(1:end-4),'_header.csv')
starttime_file_txt=strcat(file(1:end-4),'_starttime.txt')

%csvwrite(target_file_csv,frameTimestamps_sec)

fid = fopen(starttime_file_txt,'wt');
fprintf(fid, starttime);
fclose(fid);
%target_file_mat = convertStringsToChars(target_file_mat)
%save(target_file_mat,frameTimestamps_sec)

clear all