clear;
% INT549:thr=7
prompt={'tau','threshold_scaling'};
answer=inputdlg(prompt);

% set to folder with tiffs
root = 'D:\McCormick_Data_May27_2021-\2021\Jun\4004_210610_E266_meso_am_1\2P';
file='4004_210610_E266_0_00012_00001.tif';
fs = dir(fullfile(root,file ));
fname = fs(1).name;
fname = fullfile(root, fname);

% try
header = imfinfo(fname);
stack = loadFramesBuff(fname,1,1,1);

% get relevant infomation for TIFF header
artist_info = header(1).Artist;
% retrieve ScanImage ROIs information from json-encoded string
artist_info = artist_info(1:find(artist_info == '}', 1, 'last'));
artist = jsondecode(artist_info);

hSIh = header(1).Software;
hSIh = regexp(splitlines(hSIh), ' = ', 'split');
for n=1:length(hSIh)
    if strfind(hSIh{n}{1}, 'SI.hRoiManager.scanVolumeRate')
        fs = str2double(hSIh{n}{2});
    end
end

%%
si_rois = artist.RoiGroups.imagingRoiGroup.rois;

% get ROIs dimensions for each z-plane
nrois_tot = numel(si_rois);
Ly_all = [];
Lx_all = [];
cXY_all = [];
szXY_all = [];

for k = 1:nrois_tot
    zs(k,1)=si_rois(k).zs;
    Ly_all(k,1) = si_rois(k).scanfields(1).pixelResolutionXY(2);
    Lx_all(k,1) = si_rois(k).scanfields(1).pixelResolutionXY(1);
    cXY_all(k, [2 1]) = si_rois(k).scanfields(1).centerXY;
    szXY_all(k, [2 1]) = si_rois(k).scanfields(1).sizeXY;
end
nplanes = numel(unique(zs));

cXY_all = cXY_all - szXY_all/2;
cXY_all = cXY_all - min(cXY_all, [], 1);
mu = median([Ly_all, Lx_all]./szXY_all, 1);
imin_all = cXY_all .* mu;

[zs_unique,~,idxx]=unique(zs,'stable');
nrois_tot=accumarray(idxx,1);
[nrois,nrois_max_idx]=max(nrois_tot);
nrois_max_idx=find(zs==zs_unique(nrois_max_idx));
n_rows_sum = sum(Ly_all(nrois_max_idx));
n_flyback = (size(stack, 1) - n_rows_sum) / max(1, (nrois - 1));

if numel(unique(nrois_tot))==1
    f=msgbox('Length of each plane is same!!!');
else
    f=msgbox('Length of each plane is different. Seperate planes!!!');
end

for n=1:nplanes
    idx=find(idxx==n);
    Ly=Ly_all(idx);
    Lx=Lx_all(idx);
    szXY=szXY_all(idx);
    cXY=cXY_all(idx);
    imin=imin_all(idx,:);
    nrois=numel(Ly);
    
    irow = [0 cumsum(Ly'+n_flyback)];
    irow(end) = [];
    irow(2,:) = irow(1,:) + Ly';
    
    data{n}.fs = fs;
    data{n}.nplanes = 1;
    data{n}.nrois = size(irow,2);
    
    res=round(header(1).XResolution)/10^4;
    diameter = res*15;
    data{n}.diameter = round(diameter);
    if isfloat(diameter)
        f=msgbox(['Diameter ',num2str(diameter),' rounded to ',num2str(data{n}.diameter)]);
    end
    
    dx=[]; dy=[];
    data{n}.nrois = size(irow,2);
    for i = 1:size(irow,2)
        dx(i) = int32(imin(i,2));
        dy(i) = int32(imin(i,1));
        data{n}.lines{i} = irow(1,i):(irow(2,i)-1);
    end
    data{n}.dx=dx-min(dx);
    data{n}.dy=dy-min(dy);
    
    %% SAVE JSON
    tau=str2double(answer{1});
    thr=str2double(answer{2});
    data{n}.mesoscan = 1;
    data{n}.tau = tau;
    data{n}.data_path = [];
    data{n}.save_path0 = [];
    
    data{n}.mesoscan = 1;
    data{n}.delete_bin = 1;
    data{n}.save_mat=1;
    data{n}.combined = 1;
    data{n}.reg_tif=1;
    data{n}.nimg_init=300;
    data{n}.batch_size=500;
    data{n}.maxregshift=0.01;
    data{n}.keep_movie_raw = 0;
    data{n}.nonrigid = 1;
    data{n}.threshold_scaling=thr;
    data{n}.high_pass=150;
    data{n}.num_workers_roi = 5;
    
    d = jsonencode(data{n});
    fileID = fopen(fullfile(root, [file(1:end-4),'_P',num2str(n),'ops.json']),'w');
    fprintf(fileID, d);
    fclose(fileID);
end