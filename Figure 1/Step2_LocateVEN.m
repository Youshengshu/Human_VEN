Fidx=dir('MAX_NeuN_20240905_S15_section2_cVEN4.tif')
zStep=1;%um
for k=1:length(Fidx)
    filename=Fidx(k).name;    
    disp([num2str(k),'/',num2str(length(Fidx)),'--',filename])
    figure(300),clf
    imageData =imread(filename);
%     imageData=imadjust(imageData);
    imshow(imageData);hold on
    
    % 定位pia
    
    [x,y]=ginput(2);
    fit=polyfit(x,y,1);
    xx=1:1:size(imageData,2);
    plot(xx,polyval(fit,xx));
    % idx=find(coloc(:,8)==0);
    % coloc_DAPIx=coloc(idx,2);
    % coloc_DAPIy=coloc(idx,3);
    % coloc_DAPIy=coloc_DAPIy-polyval(fit,coloc_DAPIx);
 
%     disp('定位VEN坐标，回车键确认');
%     VEN_location=[];
%     VEN_location=ginput;
%     disp('VEN_location=');
%     disp(VEN_location);
%     VEN_x = VEN_location(:,1);
%     VEN_y = VEN_location(:,2);
%     VEN_y=VEN_y-polyval(fit,VEN_x);
% 
%     load([filename(10:end-4),'.mat'])
%     save([filename(10:end-4),'.mat'],'-append','VEN_location','VEN_y');
    
%     disp('定位TRI坐标，回车键确认');
%     TRI_location=[];
%     TRI_location=ginput;
%     disp('TRI_location=');
%     disp(TRI_location);
%     TRI_x = TRI_location(:,1);
%     TRI_y = TRI_location(:,2);
%     TRI_y=TRI_y-polyval(fit,TRI_x);

     disp('定位ETPC坐标，回车键确认');
    ETPC_location=[];
    ETPC_location=ginput;
    disp('ETPC_location=');
    disp(ETPC_location);
    ETPC_x = ETPC_location(:,1);
    ETPC_y = ETPC_location(:,2);
    ETPC_y=ETPC_y-polyval(fit,ETPC_x);

    % 获取图像Zstack信息
    filename=['NeuN_',filename(10:end)]
    disp([num2str(k),'/',num2str(length(Fidx)),'--',filename])
    info = imfinfo(filename);
    zstack=length(info);
   
    
    load([filename(6:end-4),'.mat'])
    save([filename(6:end-4),'.mat'],'-append','zstack','zStep','TRI_location','TRI_y','ETPC_location','ETPC_y');
    
end