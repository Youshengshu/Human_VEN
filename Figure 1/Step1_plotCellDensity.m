close all;clear all;

scale=1.1455; %0.5702594;
zPlane=[2:10];% 注意如果原始图片为1：n层，则该参数应填写2：n-1
zStep=1;%1
fillCutoff=0.7;%最后的参数由0.8改为0.7
%% load data
filename='20240815_S7_cVEN';

%% 3 steps to get cell position in image stacks

[images,centroids,areas,pixflags] = segmentObjectsByClustering(['DAPI_',filename,'.tif'],zPlane,4,zStep,1);
centroids3d = get3DCentroids(centroids,areas,pixflags,zStep,1,zPlane);
[coloc,chanIms,ims] = evalMultChannels(centroids3d.threeD,{['NeuN_',filename]},zPlane,1,fillCutoff);
figure(100)
plot(coloc(:,2),coloc(:,3),'ro');

figure(200)
plotZ=3;
idx=find(coloc(:,8)==1);%%找到NeuN+的DAPI centroid
NeuN_percent=length(idx)/size(coloc,1);
plot(coloc(idx,2),coloc(idx,3),'ro');
% plot3dMultiChannel(coloc(idx,:),chanIms{1,1},1,plotZ,2);

%test
% [images,centroids,areas,pixflags] = segmentObjectsByClustering('DAPIfile.tif',2:8,2,2,1);
% centroids3d = get3DCentroids(centroids,areas,pixflags,2,1,[2,8]);
% [coloc,chanIms,ims] = evalMultChannels(centroids3d.stereo,{'NeuNfile'},2:8,1,0.9);
% figure(100)
% plot(coloc(:,2),coloc(:,3),'ro');
% figure(200)
% plotZ=3;
% idx=find(coloc(:,8)==1);%%找到NeuN+的DAPI centroid
% plot(coloc(idx,2),coloc(idx,3),'ro');
% plot3dMultiChannel(coloc(idx,:),chanIms{1,1},1,plotZ,2);
%% what is saved in coloc
%%% Columns 1-7 are inherited from step #2.
%   one row per object, columns as follows:
%       1: object label number (matches column 5 from centroidsOUT.twoD)
%       2: X position of 3D centroid
%       3: Y position of 3D centroid
%       4: Z position of 3D centroid
%       5: object type 0-4; 0=truncated with widest plane at boundary of
%           stack; 1=bottom contained within stack; 2=top contained within
%           stack; 3=widest point is within stack and within object;
%           4=truncated by stack and has more than one area "zero crossing"
%       6: pixflags (max)
%       7: maximum area of 2D objects in group
%%% The final n+1 columns (n=number of non-DAPI channels between 1 and 3) will contain
%       8: positive (1)/negative (0) labels for each channel (per cell).(Type 0: unlabeled (DAPI only)  Type 1: Channel 1 only)
%       9: the joint labeling category (0-7, as above) in the last column.
%% find the pia
figure(300)
imageData =imread(['NeuN_',filename,'.tif'],plotZ);
imageData=imadjust(imageData);
imshow(imageData);hold on

yDist=scale*size(imageData,1);
xDist=scale*size(imageData,2);
% 定位pia
[x,y]=ginput(2);
fit=polyfit(x,y,1);
xx=1:1:size(imageData,2);
plot(xx,polyval(fit,xx));
idx=find(coloc(:,8)==1);
coloc_NeuNx=coloc(idx,2);
coloc_NeuNy=coloc(idx,3);
coloc_NeuNy=coloc_NeuNy-polyval(fit,coloc_NeuNx);
idx=find(coloc(:,8)==0);
coloc_DAPIx=coloc(idx,2);
coloc_DAPIy=coloc(idx,3);
coloc_DAPIy=coloc_DAPIy-polyval(fit,coloc_DAPIx);



%定位皮层厚度
disp('定位皮层厚度');
[x,y]=ginput(2);
fit=polyfit(x,y,1);
xx=1:1:size(imageData,2);
plot(xx,polyval(fit,xx));
depth=diff(y);

save ([filename,'.mat'],'centroids3d','coloc','coloc_NeuNy','coloc_DAPIy','depth','scale','xDist','yDist');

%% Calculation of numerical density
Binwidth=100;
BinL=0:Binwidth:4800-Binwidth;
BinR=BinL+Binwidth;
Bins=[BinL',BinR'];
cellDens_Record=[];cellNum_Record=[];
for n=1:length(Bins)
    cellNum=length(find(coloc_NeuNy>Bins(n,1)&coloc_NeuNy<Bins(n,2)));
    cellDensity=(cellNum/(Binwidth*scale*xDist*length(zPlane)))*10^9; %cells per mm3
    cellNum_Record=[cellNum_Record,cellNum];
    cellDens_Record=[cellDens_Record,cellDensity];
    
end
figure (400)
Normalized_Dist=mean(Bins,2)'/depth;
plot(cellDens_Record,Normalized_Dist,'ro')


