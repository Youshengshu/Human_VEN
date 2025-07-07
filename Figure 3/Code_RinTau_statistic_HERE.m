clear all;clc
close all;
% cd ("G:\VEN project_KW\Data analysis\Patch-seq mutimodal\Etype\matfile_AP_Para\matfile_FI")

GroupName={'cluster1','cluster2','cluster3','cluster4'};
WantState=1;%1/2：分析with/without holding状态下的数据

colormap=[0.12,0.56,1;1,0.38,0;
    0.75,0.75,0.75;1,0.87,0.68;0.75,0.75,0.75];
%蓝色[0.12,0.56,1];浅灰色[0.75,0.75,0.75];橙色[1,0.38,0];浅橙色[1,0.87,0.68];红色0.69,0.09,0.12;

%% cluster infomation for etype
Fidx={};
totaldata=readmatrix('clusteringOutput_using.csv', 'OutputType', 'string');
total_filename=totaldata(:,1);
etype=str2double(totaldata(:,3));% cluster 4 or cluster 3
for kk=1:length(unique(etype))
    idx=[];
    idx=find(ismember(etype,kk)==1);
    Fidx{kk}=total_filename(idx);

     for k=1:length(Fidx{kk})
        filename=Fidx{kk}(k);
        filename=convertStringsToChars(filename);
        
        %% loading data
        filename=['DATA_RinTau_',filename(14:end)];
        if ~exist(filename)
            continue
        end
        load(filename)
%         if sum(strcmp(filename,totalfile))~=0
%         
%         else
%             paramRecord{kk}(k,:)=nan(1,6);
%             waveformRecord{kk}(:,k)=nan(22501,1);
%             continue
%         end
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(GroupName)),'--',filename])
        %% ------------- reading data ---------------------------
%         stateIdx{1}=find(round(Vrest)<=-60);%筛选Vrest<-60mV的trial
%         stateIdx{2}=find(abs(round(Ihold))<=2);%筛选without holding的trial
        wantidx=find(StimAmp==-100);
%         wantidx=intersect(stateIdx{WantState},wantidx);         
        if ~isempty(wantidx)
            paramRecord{kk}(k,:)=[mean(Vrest(wantidx)),mean(RinRecord(wantidx)),nanmean(tauRecord(wantidx)),mean(sagRatio(wantidx)),mean(Cm(wantidx)),mean(Ihold)];
            waveformRecord{kk}(:,k)=mean(waveformV(:,wantidx)-repmat(mean(waveformV(tspanV<0,wantidx)),length(tspanV),1),2);
            figure(1),hold on
            subplot(1,length(GroupName),kk)
            tspanV_downsampled=downsample(tspanV,4);
            waveform_downsampled=downsample(waveformRecord{kk}(:,k),4);
            plot(tspanV_downsampled,waveform_downsampled,'color',colormap(kk,:))
            ylim([-10,4])
            box off
%             FilenameRecord{kk}{k,1}=filename;
%             pause
        else
            paramRecord{kk}(k,:)=nan(1,6);
            waveformRecord{kk}(:,k)=nan(45001,1);
        end
        
    end
%     pause
end
ParamName={'Vrest','Rin','tau','sagRatio','Cm'};
% save('Results_Rintau.mat','paramRecord','ParamName',...
%     'Fidx','GroupName','waveformRecord','tspanV','Fidx');

filenameRecord=Fidx{1};
for i = 1:length(filenameRecord)
    str = filenameRecord{i}; % 获取当前字符串
    extractedNames{i} = str(14:end-4); % 提取第14到倒数第4个字符
end
extractedNames=transpose(extractedNames);
RinRecord=paramRecord{1}(:,2);
passiveVmRecord=waveformRecord{1};
save('Results_Rintau.mat','extractedNames','RinRecord','passiveVmRecord');



totalpara=[];
totalfilename=[];
for kk=1:length(paramRecord)
    totalpara=[paramRecord{kk};totalpara];
    totalfilename=[Fidx{kk};totalfilename];
end
%%
% clear all;
load('Results_Rintau.mat')
figure(2),clf
for m=1:length(ParamName)
    subplot(2,7,m)
    tag=[];
    Data=[];
    for kk=1:length(Fidx)
        predata=paramRecord{kk}(:,m);
        predata(predata==0)=nan;
        data{kk}=predata;        
        N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
    
        Data=[Data;data{kk}];
        tag=[tag;kk*ones(size(data{kk}))];
          % 存储每组的均值、标准误和样本量
        group_means(kk) = Mean;
        group_ses(kk) = Sem;
        group_ns(kk) = N;
    end
        % 绘制小提琴图
    violinplot(Data, tag); hold on
    
    % 输出每组的均值、标准误和样本量
    fprintf('\nParameter: %s\n', ParamName{m});
    for kk = 1:length(Fidx)
        fprintf('Group %d: Mean = %.4f, SE = %.4f, n = %d\n', kk, group_means(kk), group_ses(kk), group_ns(kk));
    end
    
    % 组间差异检验
    [p, tbl, stats] = kruskalwallis(Data, tag, 'off');
    title([ParamName{m}, ' (p = ', num2str(round(p, 4)), ')'])
    
    % 如果整体差异显著，进行两两比较
    if p < 0.05
        [c, ~, ~, group_names] = multcompare(stats, 'Display', 'off');
        fprintf('\nMultiple comparisons for %s:\n', ParamName{m});
        for i = 1:size(c, 1)
            fprintf('Group %s vs Group %s: p = %.4f\n', ...
                group_names{c(i, 1)}, group_names{c(i, 2)}, c(i, 6));
        end
    end
    
    % 设置图形属性
    ylabel(ParamName{m})
    set(gca, 'PositionConstraint', 'innerposition')
    box off
 
end



figure(4),clf
for kk=1:length(Fidx)
    plot(tspanV,nanmean(waveformRecord{kk},2),'color',colormap(kk,:)),hold on
end


