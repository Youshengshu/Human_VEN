clear all;clc
Fidx{1}=dir('2023*.mat');
Fidx{2}=dir('2024*.mat');

Binwidth=0.03;
BinL=0:Binwidth:1.2-Binwidth;
BinR=BinL+Binwidth;
Bins=[BinL',BinR'];

for kk=1:length(Fidx)
    TRI_Record=[];
    Depth=[];
    for k=1:length(Fidx{kk})
        TRI_y=[];
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),filename])
        
        for n=2:length(Bins)
            % NeuN density
            coloc_NeuNy_scaled=coloc_NeuNy./depth;
            cellNum=length(find(coloc_NeuNy_scaled>Bins(n,1)&coloc_NeuNy_scaled<Bins(n,2)));
            cellDensity=(cellNum/(Binwidth*depth*xDist*zstack*zStep))*10^9; %cells per mm3
            cellNum_Record{1}(k,n)=cellNum;
            cellDens_Record{1}(k,n)=cellDensity;
            % DAPI density
            coloc_DAPIy_scaled=coloc_DAPIy./depth;
            cellNum=length(find(coloc_DAPIy_scaled>Bins(n,1)&coloc_DAPIy_scaled<Bins(n,2)));
            cellDensity=(cellNum/(Binwidth*depth*xDist*zstack*zStep))*10^9; %cells per mm3
            cellNum_Record{2}(k,n)=cellNum;
            cellDens_Record{2}(k,n)=cellDensity;
            
        end
        
        % VEN location
        %     VEN_y_scaled=VEN_y./depth;
        %     VEN_Record=[VEN_Record;VEN_y_scaled];
        %
        %     figure(1),clf
        %     plot(cellDens_Record{1},mean(Bins,2)','ro'),hold on
        %     plot(cellDens_Record{2},mean(Bins,2)','bo')
        
        % TRI location
        TRI_y_scaled=TRI_y./depth
        TRI_Record=[TRI_Record;TRI_y_scaled];
        
        figure(1),clf
        plot(cellDens_Record{1},mean(Bins,2)','ro'),hold on
        plot(cellDens_Record{2},mean(Bins,2)','bo')
        
        % ETPC location
        %     ETPC_y_scaled=ETPC_y./depth
        %     ETPC_Record=[ETPC_Record;ETPC_y_scaled];
        %
        %     figure(1),clf
        %     plot(cellDens_Record{1},mean(Bins,2)','ro'),hold on
        %     plot(cellDens_Record{2},mean(Bins,2)','bo')
        
        %     pause
        Depth=[Depth;depth*scale];
    end
    Depth_Record{kk}=Depth;
end


%% location of different cell types
figure(400),clf
t = tiledlayout(1,2); % 创建一个tiledlayout
% 第一个坐标系
nexttile
% colormap=[0.5,0.54,0.53;1,0.38,0];
colormap={'-bo','-ro'};
for m=1:size(cellDens_Record,2)
    Normalized_Dist=mean(Bins,2)';
    y=nanmean(cellDens_Record{m},1);
    n=size(cellDens_Record{m},2);
    sem=nanstd(cellDens_Record{m},[],1)./sqrt(n);

    % errorbar(Normalized_Dist,y,sem,'-o','color',colormap(m,:));
    shadedErrorBar(Normalized_Dist,y,sem,colormap{m},0.2)
    view(90,90)
    hold on
    box off
    
end

% plot VEN location probability
nexttile
[y,x]=ksdensity(TRI_Record,linspace(0,max(Normalized_Dist),40),'Support','positive','Bandwidth',0.1);
plot(x,y);
view(90,90)
hold on
box off

% plot VEN location
a = 1;b = 2;
r = a + (b-a) * rand(length(TRI_Record),1);
plot(TRI_Record,r,'ro');

%% comparison of cortex depth
figure(3),clf
data={};
colormap=[0.12,0.56,1;0.75,0.75,0.75
    1,0.38,0;1,0.87,0.68;0.69,0.09,0.12];
ParamName={'depth'};
for m=1:length(ParamName)
    subplot(2,1,1)
    Data=[];
    tag=[];
    for kk=1:length(Fidx)
        data{kk}=Depth_Record{kk};
        
        N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        Data=[Data;data{kk}];
        tag=[tag;kk*ones(size(data{kk}))];
    end
    boxplot(Data,tag,'widths',0.3,'Colors','k'),hold on
    plot(0.8+0.2*rand(size(data{1})),data{1},'bo','markerfacecolor',colormap(1,:),'MarkerEdgeColor','none', 'markersize',6,'lineWidth',1),hold on
    plot(1.8+0.2*rand(size(data{2})),data{2},'ro','markerfacecolor',colormap(2,:),'MarkerEdgeColor','none', 'markersize',6,'lineWidth',1),hold on
    ylabel(ParamName{m})
%     violinplot(Data,tag);hold on
%     ylabel(ParamName{m})
    set(gca,'PositionConstraint','innerposition')
    box off
    
            %%========== test============
            A=data{1};B=data{2};
            if swtest(A,0.05)==0&swtest(B,0.05)==0
                [~,p]=ttest2(A,B);
                tmethod='two samples T-test';
            else
                [p,~]=ranksum(A,B);
                tmethod='Wilcoxon rank sum ';% or Mann-Whitney U-test
            end
            %     title({['No ',num2str(nanmean(A)),'\pm',num2str(nanstd(A)./sqrt(length(A)))]
            %         ['YES ',num2str(nanmean(B)),'\pm',num2str(nanstd(B)./sqrt(length(B)))]
            %         [tmethod]
            %         ['  p = ',num2str(p)]})
            title(['  p = ',num2str(roundn(p,-4))])
    
end

