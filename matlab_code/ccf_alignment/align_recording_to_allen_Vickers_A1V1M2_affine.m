function tform = align_recording_to_allen_Vickers_A1V1M2_affine(image,alignareas,affine)
%% Align recording to Allen CCF given an image 
% INPUTSa
% image - this could be Kmeans map, or a mean of the gcamp activity
%
% alignareas - additional areas in the form 'L area', with a space between
% L/R and area name. the user should be fairly certain of identifying these areas.
%
% affine - do an affine transformation or not. default: not an affine
% transformation. The alternative is a non-reflective similarity, aka
% rotation, translation, and scaling.
%
% OUTPUT
% tform - transformation details that can be run with imwarp for subsequent frames.

% Author: Shreya Saxena (2018)

% Modified: Evan Vickers, Jan-Feb 2021
%%
if nargin<2, alignareas={}; end
if nargin<3, affine=0; end

load('allenA1V1M2DorsalMap.mat');
         
OBrefs=[190.2623  101.4249;    % Base of L OB
        396.2620  101.4249];     % Base of R OB
    
OBRSrefs=[293.9059   97.5624;   % Center base of OBs
          293.9059  431.0245];  % Base of RS

Bregref=[293.9059  244.3372]; % Bregma Center

mean_points_refs=zeros(length(alignareas),2);
for i=1:length(alignareas)
    sidestr=alignareas{i}(1);
    try
        areaid=dorsalMapsA1V1M2.labelTable{strcmp(dorsalMapsA1V1M2.labelTable.abbreviation,alignareas{i}(3:end)),'id'};
        disp(areaid);
    catch
        error('Please name additional areas in the form L/R area, with a space between L/R and area name');
    end
    dimsmap=size(dorsalMapsA1V1M2.dorsalMapScaled);
    [x,y]=find(dorsalMapsA1V1M2.dorsalMapScaled==areaid);
    if strcmp(sidestr,'R')
        x=x(y>=dimsmap(2)/2); y=y(y>=dimsmap(2)/2);
    elseif strcmp(sidestr,'L')
        x=x(y<dimsmap(2)/2); y=y(y<dimsmap(2)/2);
    end
    mean_points_refs(i,:)=[mean(y) mean(x)];
end

disp(mean_points_refs);

figure; %suptitle('Alignment to Allen CCF');
subplot(221);
imCCF=imagesc(dorsalMapsA1V1M2.dorsalMapScaled);
imCCF.AlphaData=0.9;
axis equal off;
hold on;
for p = 1:length(dorsalMapsA1V1M2.edgeOutline)
    plot(dorsalMapsA1V1M2.edgeOutline{p}(:, 2), dorsalMapsA1V1M2.edgeOutline{p}(:, 1));
end
set(gca, 'YDir', 'reverse');
plot(OBrefs(:,1),OBrefs(:,2),'-xr','linewidth',1,'markersize',3);
htext=text(50,600,'Base of L OB, R OB','color','w','fontsize',12);
title('Click on the corresponding points on your image from left to right, then press enter','fontsize',12);

handle=subplot(222);
imReal=imshow(image); axis equal off;hold on;
%imReal.AlphaData=0.8;
[x,y] = getline(handle);
OBpoints=[x y];
plot(OBpoints(:,1),OBpoints(:,2),'-xk','linewidth',2,'markersize',3);

subplot(221);
plot(OBRSrefs(:,1),OBRSrefs(:,2),'-xr','linewidth',2,'markersize',3);
delete(htext);
htext=text(50,600,'Center base of L OB, Base of RS','color','w','fontsize',12);
title('Click on the corresponding points on your image from up to down, then press enter','fontsize',12);

handle=subplot(222);
[x,y] = getline(handle);
OBRSpoints=[x y];
plot(OBRSpoints(:,1),OBRSpoints(:,2),'-xk','linewidth',2,'markersize',3);

subplot(221);
plot(Bregref(:,1),Bregref(:,2),'xr','linewidth',2,'markersize',3);
delete(htext);
htext=text(50,600,'Bregma Center','color','w','fontsize',12);
title('Click on the corresponding point on your image, then press enter','fontsize',12);

handle=subplot(222);
[x,y] = getpts(handle);
Bregpoint=[x y];
plot(Bregpoint(:,1),Bregpoint(:,2),'xk','linewidth',2,'markersize',3);

mean_points=zeros(length(alignareas),2);
for i=1:length(alignareas)
    subplot(221);
    plot(mean_points_refs(i,1),mean_points_refs(i,2),'xr','linewidth',2,'markersize',3);
    delete(htext);
    htext=text(50,600,alignareas{i},'color','w','fontsize',12);
    title('Click on the corresponding point on your image, then press enter','fontsize',12);

    handle=subplot(222);
    [x,y] = getpts(handle);
    mean_points(i,:)=[x y];
    plot(mean_points(i,1),mean_points(i,2),'xk','linewidth',2);
end
subplot(221); title('All Done');

refpoints=[OBrefs;OBRSrefs;Bregref;mean_points_refs];
points=[OBpoints;OBRSpoints;Bregpoint;mean_points];

if affine, affinestr='affine'; else, affinestr='nonreflectivesimilarity'; end
%affinestr='pwl';
%affinestr='lwm';
tform = fitgeotrans(points,refpoints,affinestr);
rotpoints = transformPointsForward(tform,points);
image(isnan(image))=0;
imagereg = imwarp(image,tform,'OutputView',imref2d(size(dorsalMapsA1V1M2.dorsalMapScaled)));

subplot(223);
plot(refpoints(:,1),refpoints(:,2),'xr',points(:,1),points(:,2),'xk',rotpoints(:,1),rotpoints(:,2),'ob','linewidth',2); 
set(gca, 'YDir', 'reverse','fontsize',12); title('Control points','fontsize',16);
legend('Allen Control Points','Unaligned Points','Points post-alignment')

subplot(224)
imReg=imshow(imagereg);
%imReg.AlphaData=0.8;
axis equal off; hold on;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMapsA1V1M2.edgeOutline{p}(:, 2), dorsalMapsA1V1M2.edgeOutline{p}(:, 1));
end
plot(refpoints(:,1),refpoints(:,2),'xr',rotpoints(:,1),rotpoints(:,2),'ob','linewidth',2); 
set(gca, 'YDir', 'reverse');

figure;imshow(imagereg)
axis equal off; hold on;
for p = 1:length(dorsalMaps.edgeOutline)
    plot(dorsalMapsA1V1M2.edgeOutline{p}(:, 2), dorsalMapsA1V1M2.edgeOutline{p}(:, 1));
end
plot(refpoints(:,1),refpoints(:,2),'xr',rotpoints(:,1),rotpoints(:,2),'ob','linewidth',2); 
set(gca, 'YDir', 'reverse');
outFile = input('Enter output file name:','s')
saveas(gcf,strcat(outFile,'_A1V1M2_CCF_affine.png'))

imCCF=imread(strcat(outFile,'_A1V1M2_CCF_affine.png'));
invtform=invert(tform);
imCCFinv=imwarp(imCCF,invtform);
figure;imshow(imCCFinv)
saveas(gcf,strcat(outFile,'_A1V1M2_CCF_affine_inverse.png'))

save(strcat(outFile,'A1V1M2_affine'))
end