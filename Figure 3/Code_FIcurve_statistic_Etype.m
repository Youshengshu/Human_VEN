clear all;clc
close all;

cd ("G:\VEN project_KW\Data analysis\Patch-seq mutimodal\Etype\matfile_AP_Para\matfile_FI")
Fidx=dir('DATA_FIcurve*.mat');
WantState=1;%1/2：分析with/without holding状态下的数据
Num=6;%选择大部分细胞都会达到的一个AP Num

%Binwidth=10;
Binwidth=20;
BinL=0:Binwidth:80;
BinR=BinL+Binwidth;
BinsFIcurveS=[BinL',BinR'];

Binwidth=50;
BinL=100:Binwidth:300;
BinR=BinL+Binwidth;
BinsFIcurveM=[BinL',BinR'];

Binwidth=200;
BinL=400:Binwidth:800;
BinR=BinL+Binwidth;
BinsFIcurveL=[BinL',BinR'];

Binwidth=500;
BinL=1000:Binwidth:4000;
BinR=BinL+Binwidth;
BinsFIcurvemax=[BinL',BinR'];

BinsFIcurve{1}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL;BinsFIcurvemax];
BinsFIcurve{2}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL;BinsFIcurvemax];
BinsFIcurve{3}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL;BinsFIcurvemax];
BinsFIcurve{4}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL;BinsFIcurvemax];
BinsFIcurve{5}=[BinsFIcurveS;BinsFIcurveM;BinsFIcurveL;BinsFIcurvemax];

%%
colormap=[0.12,0.56,1;1,0.38,0;
    0.75,0.75,0.75;1,0.87,0.68;0.75,0.75,0.75];
%蓝色[0.12,0.56,1];浅灰色[0.75,0.75,0.75];橙色[1,0.38,0];浅橙色[1,0.87,0.68];红色0.69,0.09,0.12;
%1: cVEN;2:hugePC;3:PC;4:VEN;5:unhealthy cell



%% cluster infomation for ttype
% Cluster=readtable('cluster_ttype.xlsx');
% Matfile=readtable('matfile.xlsx');
% clusterType=unique(Cluster.cluster);
%
% Fidx={};
% sampe_cluster=[];
% for kk=1:length(clusterType)
%     iidxRecord=[];
%     idx=find(Cluster.cluster==clusterType(kk));
%     sampe_cluster=Cluster.sampleName(idx);
%     for i=1:length(sampe_cluster)
%         iidx=find(strcmp(Matfile.sampleName,sampe_cluster(i)));
%         iidxRecord=[iidxRecord;iidx];
%     end
%     Fidx{kk}=cellfun(@(x)strcat('DATA_FIcurve_',x),Matfile.matfile(iidxRecord),'UniformOutput',false);
% end
%% cluster infomation for etype
Fidx={};
totaldata=readmatrix('clusteringOutput_using.csv', 'OutputType', 'string');
total_filename=totaldata(:,1);
etype=str2double(totaldata(:,2));% cluster 4 or cluster 3
for kk=1:length(unique(etype))
    idx=[];
    idx=find(ismember(etype,kk)==1);
    Fidx{kk}=total_filename(idx);
    
    for k=1:length(Fidx{kk})
        filename=Fidx{kk}(k);
        filename=convertStringsToChars(filename);
        load(filename)
        %         filename=cell2mat(Fidx{kk}(k));
        %         load([filename,'.mat'])
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        dt=tspanV(2)-tspanV(1);
        
        stateIdx{1}=find(round(Vrest)<=-60);%筛选Vrest<-60mV的trial
        stateIdx{2}=find(abs(round(Ihold))<=2);%筛选without holding的trial
        %         %% fitting ReLU
        %         Stimidx=stateIdx{WantState};
        %         uniqueStim=unique(StimAmp(Stimidx));
        %% ------------- spike train para---------------------------
        uniqueStim=unique(StimAmp);
        if ~isempty(uniqueStim)
            StimAmp_act=zeros(size(StimAmp))';
            spkFreq_act=zeros(size(StimAmp))';
            for i=1:length(uniqueStim)
                idx=find(StimAmp==uniqueStim(i));
                StimAmp_act(i)=mean(StimAmp(idx));
                spkFreq_act(i)=round(mean(spkFreq(idx)));
            end
            
            Fun=@(a,b,x) max([zeros(length(x),1),a.*x(:)+b],[],2);
            fo= fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0.001,-1000],...
                'Upper',[10,1000], ...
                'Startpoint',[0.1,-10]);
            
            [tempx,idx]=sort(StimAmp_act);
            tempy=spkFreq_act(idx);
            tempy(end-2:end)=[];tempx(end-2:end)=[];
            idx=find(tempy<=25);
            
            x=[0;tempx(idx)];
            y=[0;tempy(idx)];
            
            [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
            xx=linspace(0,x(end),100);
            yy=feval(Model,xx);
            
            Fre_threshold=max(spkFreq(find(StimAmp==min(rheobase))));
            % Fre_threshold=max(Instant_spkFreq(find(StimAmp==min(rheobase),1)));

            maxFre=max(spkFreq);
            
            spkparamRecord{kk}{11}(:,k)=-(Model.b/Model.a);
            spkparamRecord{kk}{12}(:,k)=Model.a;
            spkparamRecord{kk}{13}(:,k)=Fre_threshold;
            spkparamRecord{kk}{14}(:,k)=maxFre;
        end
        figure(1),clf
        plot(StimAmp_act,spkFreq_act,'o','markersize',6),hold on
        plot(xx,yy,'m')
        xlim([0,400])
        ylim([0,30])
        box off
        pause(0.1)
        
        %% ------------- Bin FIcurve ---------------------------
        for i=1:size(BinsFIcurve{1},1)
            idx=find(StimAmp>=BinsFIcurve{1}(i,1)&StimAmp<BinsFIcurve{1}(i,2));
            idx=intersect(stateIdx{WantState},idx);
            if ~isempty(idx)
                spkFreqBin{kk}(k,i)=mean(spkFreq(idx));
            else
                spkFreqBin{kk}(k,i)=nan;
            end
        end
        
        %% ------------- Bin spkshape,1~5 AP first spike ---------------------------
        Wantidx=find(ismember(spkFreq,[Num*2:2:Num*2+20]));
        Wantidx=intersect(stateIdx{WantState},Wantidx);
        [~,iidx]=min(StimAmp(Wantidx));
        
        if ~isempty(Wantidx)
            Wantidx=Wantidx(iidx);
            spkwaveform_temp=[];
            HWtemp=[];
            AHPtemp=[];
            Thresholdtemp=[];
            Amptemp=[];
            max_slopetemp=[];
            min_slopetemp=[];
            AIS_slope=[];
            Soma_slope=[];
            
            for jj=1%:length(spkTime{Wantidx})
                data=waveformV(:,Wantidx);
                [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,AIS_maxslope,Soma_maxslope,max_slope,min_slope...
                    HWidx,AHPidx, Thresholdidx,min_slope2idx] = AP_Statistic_new(data,spkTime{Wantidx}(jj),tspanV,1);
                
                
                spkwaveform_temp(:,jj)=spkwaveform;
                HWtemp(jj,:)=HW;
                AHPtemp(jj)=AHP;
                Thresholdtemp(jj)=Threshold;
                Amptemp(jj)=Amp;
                max_slopetemp(jj)=max_slope;
                min_slopetemp(jj)=min_slope;
                AIS_slope(jj)=AIS_maxslope;
                Soma_slope(jj)=Soma_maxslope;
                
%                 figure(100)
%                 subplot(round(length(Fidx{kk})/6)+1,6,k)
%                 plot(tspanV,data),hold on
%                 title([filename(17:end-4)], 'Interpreter','none')
%                 ylim([-80 60])
%                 xlim([-0.1,0.8])
%                 box off
                
                figure(200)
                H=[];
                %                 H(1)=subplot(3,length(spkTime{Wantidx}),jj);
                H(1)=subplot(3,1,1);
                plot(ttspanspk,spkwaveform_temp(:,jj)),hold on
                plot(ttspanspk(HWidx),spkwaveform_temp(HWidx,jj),'ro')
                plot(ttspanspk(AHPidx),spkwaveform_temp(AHPidx,jj),'ko')
                plot(ttspanspk(Thresholdidx),spkwaveform_temp(Thresholdidx,jj),'co')
                ylim([-70 40])
                
                
                H(2)=subplot(3,1,2);
                %                 H(2)=subplot(3,length(spkTime{Wantidx}),jj+length(spkTime{Wantidx}));
                slope(:,jj)=[nan;diff(spkwaveform_temp(:,jj))/(dt*1000)];
                slope2=[nan;diff(slope(:,jj))/(dt*1000)]; slope2=smooth(slope2,5);
                %                 slope3=[nan;diff(slope2(:,jj))/(dt*1000*100)];
                %                 slope3=smooth(slope3,5);
                
                plot(ttspanspk,slope2','b'),hold on
                plot(ttspanspk(min_slope2idx),slope2(min_slope2idx),'ro')
                %                 plot(ttspanspk,slope3','k'),hold on
                %
                
                H(3)=subplot(3,1,3);
                %                 H(3)=subplot(3,length(spkTime{Wantidx}),jj+length(spkTime{Wantidx})*2);
                plot(H(3),spkwaveform_temp(:,jj),slope(:,jj)),hold on
                axis([-70 50 -100 500])
                %                 linkaxes(H,'x')
                drawnow
                
            end
            
            %第一个spike waveform
            spkparamRecord{kk}{1}(:,k)=mean(HWtemp);
            spkparamRecord{kk}{2}(:,k)=mean(AHPtemp);
            spkparamRecord{kk}{3}(:,k)=mean(Thresholdtemp);
            spkparamRecord{kk}{4}(:,k)=mean(Amptemp);
            if  abs(round(Ihold(Wantidx)/10))==0
                spkparamRecord{kk}{5}(:,k)=Vrest(Wantidx);
            else
                spkparamRecord{kk}{5}(:,k)=nan;
            end
            spkparamRecord{kk}{6}(:,k)=min(rheobase);
            spkparamRecord{kk}{7}(:,k)=AIS_slope(1);
            spkparamRecord{kk}{8}(:,k)=Soma_slope(1);
            spkparamRecord{kk}{9}(:,k)=max_slopetemp(1);
            spkparamRecord{kk}{10}(:,k)=min_slopetemp(1);
            spkwaveformRecord{kk}(k,:)=spkwaveform_temp(:,1);
        else
            ttspanspk=[];
            spkparamRecord{kk}{1}(:,k)=nan;
            spkparamRecord{kk}{2}(:,k)=nan;
            spkparamRecord{kk}{3}(:,k)=nan;
            spkparamRecord{kk}{4}(:,k)=nan;
            spkparamRecord{kk}{5}(:,k)=nan;
            spkparamRecord{kk}{6}(:,k)=nan;
            spkparamRecord{kk}{7}(:,k)=nan;
            spkparamRecord{kk}{8}(:,k)=nan;
            spkparamRecord{kk}{9}(:,k)=nan;
            spkparamRecord{kk}{10}(:,k)=nan;
            %             spkwaveformRecord{kk}(:,k)=nan(1376,1);
        end
%               pause
    end
%     pause
end
ParamName={'HW', 'AHP', 'Threshold', 'Amp','Vm','Rheobase','AIS_slope','Soma_slope','max slope','min slope','Model.b/a','Model.a','Fre_threshold','maxFre'};
save('Results_spkWaveform.mat','Fidx',...
    'tspanV','spkwaveformRecord','spkparamRecord','ttspanspk','ParamName')

%% ==================== FIcurve ====================
data={};
y={};
figure(1),clf
for kk=1:length(Fidx)
    
    n=1:length(BinsFIcurve{kk});
    x=mean(BinsFIcurve{kk},2);
    data{kk}=spkFreqBin{kk}(:,n);
    data{kk}(find(data{kk}==0))=nan;
    y{kk}=nanmean(data{kk},1);
    N=sum(~isnan(data{kk}),1);
    sem=nanstd(data{kk},[],1)./sqrt(N);
    
    idx=find(~isnan(y{kk}));
    x=x(idx);
    y{kk}=y{kk}(idx);
    data{kk}=data{kk}(:,idx);
    sem=sem(idx);
    xx=linspace(x(1),x(end),100);
    
    %     %     ====== log-normal distribution fit
    %     FunLogn=@(muA,muB,stdA,stdB,A,x) A.*logncdf(x,muA,stdA).*(1-logncdf(x,muB,stdB));
    %     foLogn= fitoptions('Method','NonlinearLeastSquares','MaxIter',10000,'TolX',1e-20,'TolFun',1e-20,...
    %         'Lower',[log(1),log(200),0,0,2],...
    %         'Upper',[log(2000),log(10000),log(1000),log(1000),1000], ...
    %         'Startpoint',[log(300),log(1300),log(100),log(1000),50]);
    %
    %     [ModelLogn,~]=fit(x(:),y{kk}(:),fittype(FunLogn),foLogn)
    %     yy=feval(ModelLogn,xx);
    %% fit
    
    [ModelLogn,~]=fit(x(:),y{kk}(:),'poly2')
    yy=feval(ModelLogn,xx);
    plot(xx,yy,'-','color',colormap(kk,:));hold on
    errorbar(x,y{kk},sem,'ro','color',colormap(kk,:),'markerfacecolor',colormap(kk,:),'markersize',8)
%     pause
end
ylabel('Frequency (Hz)')
xlabel('Current injection (pA)')
box off
%
% for i=1:size(data{1},2)
%     A=data{1}(:,i);
%     B=data{2}(:,i);
%     if ~isempty(A)&~isempty(B)
%         if swtest(A,0.05)==0&swtest(B,0.05)==0
%             [~,p]=ttest2(A,B);
%         else
%             [p,~]=ranksum(A,B);% or Mann-Whitney U-test
%         end
%     end
%
%     figure(1)
%     hold on
%     if p<0.05&p>=0.01;
%         text(x(i),y{1}(i)+20,'*','fontsize',20);
%     elseif p<0.01&p>=0.001;
%         text(x(i),y{1}(i)+20,'**','fontsize',20);
%     elseif p<0.001
%         text(x(i),y{1}(i)+20,'***','fontsize',20);
%     end
% end

%% ====================ReLU fitting statistic ====================
% figure(4),clf
% figure(400),clf
% data={};
% for k=1:length(FIcurveParamName)
%     for kk=1:length(Fidx)
%         data{kk}=FIcurveParamRecord{kk}(:,k);
%         N=sum(~isnan(data{kk}));
%         Mean=nanmean(data{kk});
%         Sem=nanstd(data{kk})./sqrt(N);
%         figure(4)
%         subplot(1,length(FIcurveParamName),k)
%         bar(kk,Mean,'facecolor',colormap(kk,:));hold on
%         errorbar(kk,Mean,Sem,'color',colormap(kk,:))
%
%         if k==1
%             figure(400)
%             subplot(length(FIcurveParamName),length(Fidx),kk)
%             hist(data{kk},10)
%         end
%     end
%     ylabel(FIcurveParamName{k})
%     %========== test===========
%     A=data{1};B=data{end};
%     if swtest(A,0.05)==0&swtest(B,0.05)==0
%         [~,p]=ttest2(A,B);
%         tmethod='two samples T-test';
%     else
%         [p,~]=ranksum(A,B);
%         tmethod='Wilcoxon rank sum ';% or Mann-Whitney U-test
%     end
%     title({['No ',num2str(nanmean(A)),'\pm',num2str(nanstd(A)./sqrt(length(A)))]
%         ['YES ',num2str(nanmean(B)),'\pm',num2str(nanstd(B)./sqrt(length(B)))]
%         [tmethod]
%         ['  p = ',num2str(p)]})
%
% end

%% ========AP waveform parameter==========
figure(5), clf
data = {};
stats_result = struct();  % 初始化统计结果结构体

for m = 1:length(ParamName)
    subplot(2, length(ParamName) / 2, m)
    Data = [];
    tag = [];
    group_means = [];
    group_ses = [];
    group_ns = [];

    for kk = 2:3%1:length(Fidx)
        data{kk} = spkparamRecord{kk}{m}';
        N = sum(~isnan(data{kk}));
        Mean = nanmean(data{kk});
        Sem = nanstd(data{kk}) ./ sqrt(N);

        Data = [Data; data{kk}];
        tag = [tag; kk * ones(size(data{kk}))];

        % 存储每组的均值、标准误和样本量
        group_means(kk) = Mean;
        group_ses(kk) = Sem;
        group_ns(kk) = N;
    end

    % 绘制小提琴图
    violinplot(Data, tag); hold on
    box off

    % 输出每组统计量
    fprintf('\nParameter: %s\n', ParamName{m});
    for kk = 2:3%1:length(Fidx)
        fprintf('Group %d: Mean = %.4f, SE = %.4f, n = %d\n', kk, group_means(kk), group_ses(kk), group_ns(kk));
    end

    % % Kruskal-Wallis 检验 (非参数检验方法，用于比较三个及以上独立样本)
    % [p, tbl, stats] = kruskalwallis(Data, tag, 'off');
    % title([ParamName{m}, ' (p = ', num2str(round(p, 4)), ')'])

    % 提取两组数据（比较两类VEN的差异）
    group1 = Data(tag == 2);
    group2 = Data(tag == 3);

    % 检验方差齐性
    [h_var, p_var] = vartest2(group1, group2);

    % 根据方差检验结果选择检验方法
    if h_var == 0  % 方差齐，使用双样本 t 检验
        [~, p] = ttest2(group1, group2);
        test_name = 't-test';
    else           % 方差不齐，使用非参数秩和检验
        [p, ~, ~] = ranksum(group1, group2);
        test_name = 'ranksum';
    end

    % 绘图标题显示参数名和 p 值
    title([ParamName{m}, ' (', test_name, ', p = ', num2str(round(p, 4)), ')'])



    % % 保存统计数据到结构体
    stats_result(m).ParamName = ParamName{m};
    stats_result(m).group_means = group_means;
    stats_result(m).group_ses = group_ses;
    stats_result(m).group_ns = group_ns;
    stats_result(m).p = p;
    % 
    % % 多重比较
    % if p < 0.05
    %     [c, ~, ~, group_names] = multcompare(stats, 'Display', 'off');
    %     comparisons = cell(size(c, 1), 1);
    %     for i = 1:size(c, 1)
    %         comparisons{i} = struct('group1', group_names{c(i, 1)}, ...
    %                                 'group2', group_names{c(i, 2)}, ...
    %                                 'p', c(i, 6));
    %     end
    %     stats_result(m).multiple_comparisons = comparisons;
    % else
    %     stats_result(m).multiple_comparisons = {};
    % end
end

%% 👇 写入 CSV 文件（确保结构体已填满）
% 第一张表：参数统计总览
summary_table = table();
for m = 1:length(stats_result)
    for kk = 1:length(stats_result(m).group_means)
        summary_table = [summary_table; {
            stats_result(m).ParamName, ...
            kk, ...
            stats_result(m).group_means(kk), ...
            stats_result(m).group_ses(kk), ...
            stats_result(m).group_ns(kk), ...
            stats_result(m).p
        }];
    end
end
summary_table.Properties.VariableNames = {'Parameter', 'Group', 'Mean', 'SE', 'N', 'p'};
writetable(summary_table, 'VENsubtypes_comparision.csv');

% 第二张表：多重比较
comp_table = table();
for m = 1:length(stats_result)
    comps = stats_result(m).multiple_comparisons;
    for i = 1:length(comps)
        comp_table = [comp_table; {
            stats_result(m).ParamName, ...
            comps{i}.group1, ...
            comps{i}.group2, ...
            comps{i}.p
        }];
    end
end
if ~isempty(comp_table)
    comp_table.Properties.VariableNames = {'Parameter', 'Group1', 'Group2', 'p_value'};
    writetable(comp_table, 'multiple_comparisons.csv');
end



%% ====================  Waveform display ====================
xlim_max=0.008;
xlim_min=-0.0015;
figure(7),clf

for kk=1:length(unique(etype))
    tempidx=find(ttspanspk>-0.004&ttspanspk<0.0040);
    for j=1:size(spkwaveformRecord{kk},1)
        data=spkwaveformRecord{kk}(j,:);
        meandata=squeeze(nanmean(spkwaveformRecord{kk},1));
        dt=ttspanspk(2)-ttspanspk(1);
%         h1=subplot(2,2,1);
%         plot(data(tempidx(2:end)),diff(data(tempidx))./(dt*1000),'--','color',colormap(kk,:)),hold on
%         xlabel('Vm (mV)');ylabel('Slope');
%         box off
        h2=subplot(2,2,2);
        plot(meandata(tempidx(2:end)),diff(meandata(tempidx))./(dt*1000),'--','color',colormap(kk,:)),hold on
        xlabel('Vm (mV)');ylabel('Slope');
        box off
%         h3=subplot(2,2,3);
%         plot(ttspanspk,spkwaveformRecord{kk},'color',colormap(kk,:)),hold on
%         axis([xlim_min,xlim_max,-65,60])
%         xlabel('Time (s)');ylabel('Vm (mV)');
%         box off
        h4=subplot(2,2,4);
        plot(ttspanspk,meandata,'color',colormap(kk,:)),hold on
        axis([xlim_min,xlim_max,-65,60])
        xlabel('Time (s)');ylabel('Vm (mV)');
        box off
    end
end
