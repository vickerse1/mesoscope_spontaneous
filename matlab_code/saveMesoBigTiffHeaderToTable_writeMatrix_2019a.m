prompt='Enter path: '
path=input(prompt,'s');
path=strcat(path,'\')

prompt='Enter filename: '
filename=input(prompt,'s');

file=strcat(path,filename);

[header,Aout,imgInfo,rawStream] = scanimage.util.opentif(file);

frameTimestamps_sec=header.frameTimestamps_sec;

target_file=strcat(file(1:end-4),'_header.csv')

writematrix(frameTimestamps_sec,target_file);

%clear all