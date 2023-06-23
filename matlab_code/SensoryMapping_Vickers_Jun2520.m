% Sensory Mapping Code
close all
clear
%% Load In Widefield Video
%Inputs
prompt ={'Enter number of WF frames collected.','Enter X size.','Enter Y size.','Enter sampling frequency.','Enter number of skipped frames.','Enter fraction to downsample pixels.'};
t = 'Input';
dims = [1 50];
definput = {'15000','640','640','50','0','0.5'};
answer = inputdlg(prompt,t,dims,definput);
TotalFrames = str2num(answer{1,1});
XSize = str2num(answer{2,1});
YSize = str2num(answer{3,1});
Fs = str2num(answer{4,1});
NSkipFrames = str2num(answer{5,1});
DSFrac = str2num(answer{6,1});

%Read in .tsm files
[FileName, folder] = uigetfile('*.tsm','Select WF .tsm file.');  % Opens UI to find .tsm file
videoWF = readTSMFile(folder,FileName, XSize,YSize,TotalFrames,NSkipFrames,DSFrac);
[dF,UdF] = Calc_dFF(videoWF,Fs); % Global 10 percentile baseline
%[ROITraces,AvgROITrace] = AvgROICC(dF);
meandf=mean(dF,3); %creates mean of entire trace
figure(3)
h=heatmap(meandf);


%% Read in Spike Files, 
[FileName, folder] = uigetfile('*.mat','Select Spike2 .mat file.');  
load(fullfile(folder,FileName))

%create timeseries variables for behavior
walk.values=(walk.values-3)*10;
ts_encoder = timeseries(walk.values,(linspace(0,(walk.length*walk.interval),walk.length))');
%ts_stim = timeseries(STIM.values,(linspace(0,(STIM.length*STIM.interval),STIM.length))');
ts_whisk = timeseries(whisk.values,(linspace(0,(whisk.length*whisk.interval),whisk.length))');

%read in posthoc pupil if you want
prompt = {'Would you like to read in posthoc pupil data? (Y/N)'};
t = 'Input';
dims = [1 25];
definput = {'Y'};
PHPanswer = inputdlg(prompt,t,dims,definput);
if PHPanswer{1,1} == 'Y' 
[FileNamePHP, folderPHP] = uigetfile('*smoothed');  % Opens UI to find Spike2 file
PH_pupil = (importdata(strcat(folderPHP,FileNamePHP)))';
PH_pupil=movmean(PH_pupil,5);
    if length(camera.times)>length(PH_pupil)
    PH_pupil(length(PH_pupil)+length(camera.times)-length(PH_pupil),:)=NaN;
    end
    ts_pupil = timeseries(PH_pupil, camera.times); %match up to pupilcam frame clock
elseif PHPanswer{1,1} ~= 'Y' 
    ts_pupil = timeseries(Pupil.values,(linspace(0,(Pupil.length*Pupil.interval),Pupil.length))');
    disp('No posthoc pupil fitting. Real-time Pupil used')
end

%downsample everything to 50Hz 
T = max(ts_encoder.time);
Resamp_r = 50; %in Hz
Resamp_p = 1/Resamp_r; %in seconds
times_resamp = 0:Resamp_p:T;

tsresam_pupil = resample(ts_pupil,times_resamp);
tsresam_encoder = resample(ts_encoder,times_resamp);
%tsresam_stim = resample(ts_stim,times_resamp);
tsresam_whisk = resample(ts_whisk,times_resamp);

%Convert pupil trace to % of Maximum pupil. 
figure ('Name','Pupil Raw Trace')
plot(ts_pupil.Data)
title('First choose maximum, then choose minimum pupil values.')
[~,y1] =(ginput(1));
[~,y2] = (ginput(1));
tsresam_pupil.Data(tsresam_pupil.Data < y2) = NaN;
tsresam_pupil.Data(tsresam_pupil.Data > y1) = NaN;

%it will now prompt you to either manually or automatically determine
%maximum pupil values. You may want to choose to manually input the max
%pupil value if the mouse did not run during this chunk of file and you
%know what the maximum/running pupil value is from the spike2 file.  

prompt = {'Would you like to automatically determine (A) or manually input (M) the max pupil value?'};
t = 'Input';
dims = [1 60];
definput = {'A'};
auto_manual_pupil = inputdlg(prompt,t,dims,definput);

if auto_manual_pupil{1,1} == 'A'

    maxpup = max(tsresam_pupil.Data);
    tsresam_pupil.Data =(tsresam_pupil.Data/maxpup)*100;
    tsresam_pupil.Data=movmedian(tsresam_pupil.Data,2000);
    figure('Name','Corrected Pupil Trace')
    plot(tsresam_pupil.Data)
    histogram(tsresam_pupil.Data)
    disp('Automatic Max Pupil Determined')
    
elseif auto_manual_pupil{1,1} == 'M'
    prompt = {'Input Max Pupil Value'};
    t = 'Input';
    dims = [1 25];
    definput = {'0.6'};
    maxpup = inputdlg(prompt,t,dims,definput);
    maxpup=str2num(cell2mat(maxpup));
    tsresam_pupil.Data =(tsresam_pupil.Data/maxpup)*100;
    tsresam_pupil.Data=movmedian(tsresam_pupil.Data,2000);
    figure('Name','Corrected Pupil Trace')
    plot(tsresam_pupil.Data)
    histogram(tsresam_pupil.Data)
    disp('Max Pupil Manually Determined')
end





%% Chop Files to Widefield Trigger

WFTime = WFframec.times;
%might need to switch CamOn and CamOff
[CamOn, CamOff, CamTime] = FindStimTimes(wfcam_tr,'WindowDiff',10e10,'CamTrig',1);
[CamWFOnset, CamWFOffsets] = CompareTimes(WFTime, CamTime(CamOn),CamTime(CamOff));
WFTime = WFTime(CamWFOnset:CamWFOnset+TotalFrames-1);
% 
% STIM=OldstmIN; %VNS
% for y = 1:length(STIM.values)
%     if STIM.values(y) > 0.1
%         STIM.values(y) = 5;
%     else
%         STIM.values(y) = 0;
%     end
% end

prompt ={'What kind of Sensory Mapping? Visual/Auditory/Whisk (V/A/W)'};
t = 'Input';
dims = [1 50];
definput = {'V'};
answer = inputdlg(prompt,t,dims,definput);


if answer{1,1} == 'V'
    for y = 1:length(photoDio.values) %Visual
        if photoDio.values(y) > 0.11
            photoDio.values(y) = 5;
        else
            photoDio.values(y) = 0;
        end
    end
    
    ts_sensorystim = timeseries(photoDio.values,(linspace(0,(photoDio.length*photoDio.interval),photoDio.length))');
    
elseif answer{1,1} == 'A'
    
    for y = 1:length(sound_ef.values) %Auditory
        if sound_ef.values(y) > 1
            sound_ef.values(y) = 5;
        else
            sound_ef.values(y) = 0;
        end
    end
    
    ts_sensorystim =  timeseries(sound_ef.values,(linspace(0,(sound_ef.length*sound_ef.interval),sound_ef.length))');
    
elseif answer{1,1} == 'W'
    
    for y = 1:length(Piezo_Co.values) %Whisk
        if Piezo_Co.values(y) > 4.93
            Piezo_Co.values(y) = 5;
        else
            Piezo_Co.values(y) = 0;
        end
    end
    
    ts_sensorystim = timeseries(Piezo_Co.values,(linspace(0,(Piezo_Co.length*Piezo_Co.interval),Piezo_Co.length))');
  
else
    disp('Invalid Inputs')
    
end

tsresam_sensorystim = resample(ts_sensorystim,times_resamp);

[stimWFOnset, stimWFOffsets] = CompareTimes(tsresam_sensorystim.time, WFTime(1),WFTime(TotalFrames));
cuttsresam_sensorystim = timeseries(tsresam_sensorystim.Data(stimWFOnset:stimWFOffsets),tsresam_sensorystim.time(stimWFOnset:stimWFOffsets));

[pupilWFOnset, pupilWFOffsets] = CompareTimes(tsresam_pupil.time, WFTime(1),WFTime(TotalFrames));
cuttsresam_pupil = timeseries(tsresam_pupil.Data(pupilWFOnset:pupilWFOffsets),tsresam_pupil.time(pupilWFOnset:pupilWFOffsets));

[whiskWFOnset, whiskWFOffsets] = CompareTimes(tsresam_whisk.time, WFTime(1),WFTime(TotalFrames));
cuttsresam_whisk = timeseries(tsresam_whisk.Data(whiskWFOnset:whiskWFOffsets),tsresam_whisk.time(whiskWFOnset:whiskWFOffsets));

[encoderWFOnset, encoderWFOffsets] = CompareTimes(tsresam_encoder.time, WFTime(1),WFTime(TotalFrames));
cuttsresam_encoder = timeseries(tsresam_encoder.Data(encoderWFOnset:encoderWFOffsets),tsresam_encoder.time(encoderWFOnset:encoderWFOffsets));

%Find Stim Onsets
if answer{1,1} == 'V'
    [pks,locs] = findpeaks(cuttsresam_sensorystim.data,'MinPeakHeight',0.5,'MinPeakDistance',Resamp_r*6); %CHANGED 20 to 10 for auditory - change back
    figure ('Name','Stim locations');
    plot(cuttsresam_sensorystim,cuttsresam_sensorystim.time(locs),pks,'or')
    title('Visual Stim Peaks')
    
elseif answer{1,1} == 'A'
    [pks,locs] = findpeaks(cuttsresam_sensorystim.data,'MinPeakHeight',0.5,'MinPeakDistance',Resamp_r*6); %CHANGED 20 to 10 for auditory - change back
    figure ('Name','Stim locations');
    plot(cuttsresam_sensorystim,cuttsresam_sensorystim.time(locs),pks,'or')
    title('Auditory Stim Peaks')
    
    elseif answer{1,1} == 'W'
    [pks,locs] = findpeaks(cuttsresam_sensorystim.data,'MinPeakHeight',0.5,'MinPeakDistance',Resamp_r*6); %Need to test this for labview/spike2 triggered whisk
    figure ('Name','Stim locations');
    plot(cuttsresam_sensorystim,cuttsresam_sensorystim.time(locs),pks,'or')
    title('Whisker Stim Peaks')
end


% %% Separate out traces based on walking and not walking 
% %pick a point where the animal isn't walking and calculate 4x the stdev of
% %that section of encoder trace
% figure ('Name','Walking Trace','units','normalized','outerposition',[0 0 1 1])
% plot(tsresam_encoder.Data) 
% title('Choose stationary segment.')
% [x1,~] =(ginput(1)); 
% [x2,~] = (ginput(1)); 
% close('Walking Trace');
% x1=round(x1); 
% x2=round(x2);
% mean_walk = mean(tsresam_encoder.data(x1:x2));
% sd_walk = std(tsresam_encoder.data(x1:x2));
% sd4_walk=(sd_walk*4);
% mean4sd=mean_walk+sd4_walk;
% %find the average walking speed BEFORE each trial
% avspeed = [];
% prewalk=Resamp_r*1;
% 
% for s = 1:length(pks)
%      avspeed(:,s) = mean(tsresam_encoder.data(locs(s)-prewalk:locs(s)));
% end
% 
% 
% 
% for i = 1:length(avspeed)
%     if avspeed(i) > mean4sd
%         this_WALK = (tsresam_encoder.Data((locs(i)-pre):(locs(i)+post)));
%         this_PUP = (tsresam_pupil.Data((locs(i)-pre):(locs(i)+post)));
%         this_WHISK = (tsresam_whisk.Data((locs(i)-pre):(locs(i)+post)));
%         this_LFP = (tsresam_lfp.Data((locs(i)-pre):(locs(i)+post)));
%         this_SU = (tsresam_SU.Data((locs(i)-pre):(locs(i)+post)));
%         WALK_Sig = [WALK_Sig, this_WALK]; 
%         PUP_Sig = [PUP_Sig, this_PUP];
%         WHISK_Sig = [WHISK_Sig, this_WHISK];
%         LFP_Sig = [LFP_Sig, this_LFP]; 
%         SU_Sig = [SU_Sig, this_SU];
%         
%     elseif  avspeed(i) <= mean4sd
%         this_WALK = (tsresam_encoder.Data((locs(i)-pre):(locs(i)+post)));
%         this_PUP = (tsresam_pupil.Data((locs(i)-pre):(locs(i)+post)));
%         this_WHISK = (tsresam_whisk.Data((locs(i)-pre):(locs(i)+post)));
%         this_LFP = (tsresam_lfp.Data((locs(i)-pre):(locs(i)+post)));
%         this_SU = (tsresam_SU.Data((locs(i)-pre):(locs(i)+post)));
%         WALK_NonSig = [WALK_NonSig, this_WALK];
%         PUP_NonSig = [PUP_NonSig, this_PUP];
%         WHISK_NonSig = [WHISK_NonSig, this_WHISK];
%         LFP_NonSig = [LFP_NonSig, this_LFP]; 
%         SU_NonSig = [SU_NonSig, this_SU];
%     else
%         disp('Error Generated. Please Check inputs/outputs to FOR Loop.') 
%     end
% end


%% Determine Time to cut files
pre=Fs*1; %1 seconds
post=Fs*1; %1 seconds

%preallocate all variables for speed
this_presensorystim=[];
this_sensorystim=[];
all_presensorystim=[];
all_sensorystim=[];

for i = 1:length(pks)
    this_presensorystim = dF(:,:,locs(i)-pre:locs(i)-1);
    all_presensorystim = cat(3, all_presensorystim,this_presensorystim);
    this_sensorystim = dF(:,:,locs(i):locs(i)+post-1);
    all_sensorystim = cat(3, all_sensorystim,this_sensorystim);  
end

mean_presensory=mean(all_presensorystim,3);
mean_sensory=mean(all_sensorystim,3);
%figure(2)
%presensory_hm=heatmap(mean_presensory,'GridVisible','off','ColorMap',jet);
%figure(3)
%sensory_hsensom=heatmap(mean_sensory,'GridVisible','off','ColorMap',jet);
mean_sensory_sub=mean_sensory-mean_presensory
figure(9)
sensory_sub_hm=heatmap(mean_sensory_sub,'GridVisible','off','ColorMap',jet);

mean_sensory_sub_rot = imrotate(mean_sensory_sub,270)
mean_sensory_sub_rot_flip = flipdim(mean_sensory_sub_rot ,2)
figure(10)
sensory_sub_rot_flip_hm=heatmap(mean_sensory_sub_rot_flip,'GridVisible','off','ColorMap',jet);
FileOut = FileName(1:end-4)
saveas(gcf,FileOut+"_df_hm.png")

figure(12)
df_gray=mat2gray(mean_sensory_sub_rot_flip)
df_gray=imadjust(df_gray)
imshow(df_gray)
saveas(gcf,FileOut+"_df_gs.png")
imwrite(df_gray,FileOut+"_df_gs.tif")
save(FileOut+"_df.mat","mean_sensory_sub_rot_flip","df_gray")

%%Add in walk filter








