
clear all;clc
close all;
Fidx{1}=dir('DATA_FIcurve*ETPC*.mat');
Fidx{2}=dir('DATA_FIcurve*cVEN*.mat');

WantState=1;  % 1/2: analyze with/without holding state data
Num=6;  % Choose a number of AP spikes most cells reach

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

% Num=6;%选择大部分细胞都会达到的一个AP Num
colormap=[0.12,0.56,1;0.75,0.75,0.75
    1,0.38,0;1,0.87,0.68;0.69,0.09,0.12];
%蓝色[0.12,0.56,1];浅灰色[0.75,0.75,0.75];橙色[1,0.38,0];浅橙色[1,0.87,0.68];红色0.69,0.09,0.12;
%1: cVEN;2:hugePC;3:PC;4:VEN;5:unhealthy cell


%% cluster infomation for etype

for kk=1:length(Fidx)
    DistRecord=[];
    for k=1:length(Fidx{kk})
        %% loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        dt=tspanV(2)-tspanV(1);

        stateIdx{1}=find(round(Vrest)<=-60);
        stateIdx{2}=find(abs(round(Ihold))<=2);

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
            fo=fitoptions('Method','NonlinearLeastSquares','Lower',[0.001,-1000],'Upper',[10,1000],...
                'Startpoint',[0.1,-10]);

            [tempx,idx]=sort(StimAmp_act);
            tempy=spkFreq_act(idx);
            tempy(end-2:end)=[]; tempx(end-2:end)=[];
            idx=find(tempy<=25);

            x=[0;tempx(idx)];
            y=[0;tempy(idx)];

            [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
            xx=linspace(0,x(end),100);
            yy=feval(Model,xx);

            Fre_threshold=max(spkFreq(find(StimAmp==min(rheobase))));
            maxFre=max(spkFreq);
          

            % spkparamRecord{kk}{12}(:,k)=-(Model.b/Model.a);
            spkparamRecord{kk}{13}(:,k)=Model.a;
            spkparamRecord{kk}{14}(:,k)=Fre_threshold;
            spkparamRecord{kk}{15}(:,k)=maxFre;
        end

        figure(1), clf
        plot(StimAmp_act, spkFreq_act, 'o', 'markersize', 6), hold on
        plot(xx, yy, 'm')
        xlim([0, 400])
        ylim([0, 30])
        box off
        pause(0.1)

        % Extract spike waveform
        for i=1:size(BinsFIcurve{1},1)
            idx=find(StimAmp>=BinsFIcurve{1}(i,1) & StimAmp<BinsFIcurve{1}(i,2));
            idx=intersect(stateIdx{WantState},idx);
            if ~isempty(idx)
                spkFreqBin{kk}(k,i)=mean(spkFreq(idx));
            else
                spkFreqBin{kk}(k,i)=nan;
            end
        end

        % Extract first AP shape parameters
        Wantidx=find(ismember(spkFreq, [Num*2:2:Num*2+20]));
        Wantidx=intersect(stateIdx{WantState}, Wantidx);
        [~,iidx]=min(StimAmp(Wantidx));

        if ~isempty(Wantidx)
            Wantidx=Wantidx(iidx);
            spkwaveform_temp=[]; HWtemp=[]; AHPtemp=[]; Thresholdtemp=[]; Amptemp=[]; max_slopetemp=[]; min_slopetemp=[];
            for jj=1%:length(spkTime{Wantidx})
                data=waveformV(:, Wantidx);
                [spkwaveform, ttspanspk, HW, AHP, Threshold, Amp, inflection_idx, initialSlope, max_slope, min_slope, ...
                    HWidx, AHPidx, Thresholdidx, min_slope2idx] = AP_Statistic_new(data, spkTime{Wantidx}(jj), tspanV, 1);

                spkwaveform_temp(:,jj)=spkwaveform;
                HWtemp(jj,:)=HW;
                AHPtemp(jj)=AHP;
                Thresholdtemp(jj)=Threshold;
                Amptemp(jj)=Amp;
                max_slopetemp(jj)=max_slope;
                min_slopetemp(jj)=min_slope;

                % Inside the loop for plotting AP shape parameters
                figure(200),clf
                H=[];

                H(1)=subplot(4, 1, 1);
                plot(ttspanspk, spkwaveform_temp(:,jj)), hold on
                plot(ttspanspk(HWidx), spkwaveform_temp(HWidx,jj), 'ro')
                plot(ttspanspk(AHPidx), spkwaveform_temp(AHPidx,jj), 'ko')
                plot(ttspanspk(Thresholdidx), spkwaveform_temp(Thresholdidx,jj), 'co')
                ylim([-70 40])

                H(2)=subplot(4, 1, 2);
                slope(:,jj) = [nan; diff(spkwaveform_temp(:,jj)) / (dt*1000)];
                plot(ttspanspk, slope(:,jj), 'b', 'LineWidth', 1.5); 
                hold on
                plot(ttspanspk(inflection_idx), initialSlope, 'ko', 'MarkerFaceColor', 'k');
                xlabel('Time (ms)')
                ylabel('dV/dt (V/s)')
                title('Slope with Inflection Point')
                box off



                % New subplot for slope2, inflection_idx, and initialSlope

                H(3) = subplot(4, 1, 3);  % Adding new subplot
                slope2=[nan; diff(slope(:,jj))/(dt*1000)]; slope2=smooth(slope2,5);
                plot(ttspanspk, slope2, 'b', 'LineWidth', 1.5);  % Plot slope2
                hold on;
                plot(ttspanspk(inflection_idx), slope2(inflection_idx), 'ro', 'MarkerFaceColor', 'r');  % Highlight inflection points

                xlabel('Time (ms)')
                ylabel('Slope (V/ms)')
                title('Slope2 with Inflection Points ')
                legend({'Slope2', 'Inflection Points'}, 'Location', 'Best')
                box off

                H(4)=subplot(4, 1, 4);
                plot(spkwaveform_temp(:, jj), slope(:, jj)), hold on
                axis([-70 50 -100 500])
                
                fitidx=Thresholdidx-6:Thresholdidx;
                slopeInitial= polyfit(spkwaveform_temp(fitidx),slope(fitidx),1);
                Fit=polyval(slopeInitial,spkwaveform_temp(fitidx));
                OnsetSlope=slopeInitial(1);
                plot(H(4),spkwaveform_temp(:,jj),slope(:,jj)),hold on
                plot(H(4),spkwaveform_temp(fitidx),Fit,'r')

            end

            % Store extracted parameters
            spkparamRecord{kk}{1}(:,k)=mean(HWtemp);
            spkparamRecord{kk}{2}(:,k)=mean(AHPtemp);
            spkparamRecord{kk}{3}(:,k)=mean(Thresholdtemp);
            spkparamRecord{kk}{4}(:,k)=mean(Amptemp);
            spkparamRecord{kk}{5}(:,k)=Vrest(Wantidx);
            spkparamRecord{kk}{6}(:,k)=min(rheobase);
            spkparamRecord{kk}{7}(:,k)=max_slopetemp(1);
            spkparamRecord{kk}{8}(:,k)=min_slopetemp(1);
            spkparamRecord{kk}{11}(:,k)=Dist;
            spkparamRecord{kk}{12}(:,k)=OnsetSlope;
            spkwaveformRecord{kk}(k,:)=spkwaveform_temp(:,1);
            if ~isnan(inflection_idx)
                spkparamRecord{kk}{9}(:,k)=initialSlope;
                spkparamRecord{kk}{10}(:,k)=ttspanspk(inflection_idx);
            else
                spkparamRecord{kk}{9}(:,k)=nan;
                spkparamRecord{kk}{10}(:,k)=nan;
            end


        else
            ttspanspk=[]; % Empty if no spikes
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
            spkparamRecord{kk}{11}(:,k)=nan;
            spkparamRecord{kk}{12}(:,k)=nan;
        end
        % pause()
         Dist=[];
    end
    Axon_Dist{kk}=DistRecord;
end

ParamName={'HW', 'AHP', 'Threshold', 'Amp','Vm','Rheobase','max_slope','min_slope','initialSlope','inflection_idx','Dist','OnsetSlope','Model.a','Fre_threshold','maxFre'};
save('Results_spkWaveform.mat', 'Fidx', 'tspanV', 'spkwaveformRecord', 'spkparamRecord', 'ttspanspk', 'ParamName')

%% inflection data exclude
% 读取 CSV 文件
filename = 'withInflection.csv';  % 替换成你的实际文件名
opts = detectImportOptions(filename);
T = readtable(filename, opts);

% 获取编号列表
ETPC_keep_idx = T.ETPC(~isnan(T.ETPC));  % 去掉空值
VEN_keep_idx  = T.VEN(~isnan(T.VEN));

% 获取原始数据（假设已存在）
ETPC_all = spkparamRecord{1}{9};  % 初始数据
VEN_all  = spkparamRecord{2}{9};

% 将不在指定编号中的元素设为 NaN
ETPC_mask = true(size(ETPC_all));
ETPC_mask(ETPC_keep_idx) = false;
ETPC_all(ETPC_mask) = NaN;

VEN_mask = true(size(VEN_all));
VEN_mask(VEN_keep_idx) = false;
VEN_all(VEN_mask) = NaN;

% 回写结果
spkparamRecord{1}{9} = ETPC_all;
spkparamRecord{2}{9} = VEN_all;



%% ========AP waveform parameter==========
figure(5), clf
data = {};
for m = 1:length(ParamName)
    subplot(2, round(length(ParamName) / 2), m)
    Data = [];
    tag = [];
    group_means = [];
    group_ses = [];
    group_ns = [];

    for kk = 1:length(Fidx)
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

%% 变为箱型图
colors = [0.2 0.6 1.0; 1.0 0.4 0.1]; % 蓝色、橙色

figure(5); clf
data = {};
for m = 1:length(ParamName)
    subplot(2, round(length(ParamName) / 2), m)
    Data = [];
    tag = [];
    group_means = [];
    group_ses = [];
    group_ns = [];

    hold on
    for kk = 1:length(Fidx)
        data{kk} = spkparamRecord{kk}{m}';
        N = sum(~isnan(data{kk}));
        Mean = nanmean(data{kk});
        Sem = nanstd(data{kk}) ./ sqrt(N);

        Data = [Data; data{kk}];
        tag = [tag; kk * ones(size(data{kk}))];

        group_means(kk) = Mean;
        group_ses(kk) = Sem;
        group_ns(kk) = N;
    end

    % 使用 boxplot 画箱型图
    boxplot(Data, tag, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5)
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), 'w', 'FaceAlpha', 0, 'EdgeColor', 'k');
    end

    % 添加彩色散点
    for kk = 1:length(Fidx)
        jitterX = kk + 0.08 * randn(size(data{kk}));
        scatter(jitterX, data{kk}, 36, colors(kk, :), 'filled', ...
                'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'k');
    end
    % 设置图形属性
    ylabel(ParamName{m})
    set(gca, 'PositionConstraint', 'innerposition')
    box off

    % 自动优化 y 轴范围
    y_min = min(Data);
    y_max = max(Data);
    y_range = y_max - y_min;
    padding = 0.1 * y_range;
    ylim([y_min - padding, y_max + padding]);

    % 显示统计值
    fprintf('\nParameter: %s\n', ParamName{m});
    for kk = 1:length(Fidx)
        fprintf('Group %d: Mean = %.4f, SE = %.4f, n = %d\n', ...
                kk, group_means(kk), group_ses(kk), group_ns(kk));
    end

    % Kruskal-Wallis检验
    [p, tbl, stats] = kruskalwallis(Data, tag, 'off');
    title([ParamName{m}, ' (p = ', num2str(round(p, 4)), ')'])

    if p < 0.05
        [c, ~, ~, group_names] = multcompare(stats, 'Display', 'off');
        fprintf('\nMultiple comparisons for %s:\n', ParamName{m});
        for i = 1:size(c, 1)
            fprintf('Group %s vs Group %s: p = %.4f\n', ...
                group_names{c(i, 1)}, group_names{c(i, 2)}, c(i, 6));
        end
    end

    ylabel(ParamName{m})
    set(gca, 'PositionConstraint', 'innerposition')
    box off
end

%% para correlation vs axon origion
figure(6),clf
for m = 1:length(ParamName)
    subplot(2, round(length(ParamName)/2), m)
    hold on

    all_x = [];
    all_y = [];
    all_c = [];

    % 收集所有组的数据（用于统一拟合）
    for kk = 1:length(Fidx)
        x = spkparamRecord{kk}{11};
        y = spkparamRecord{kk}{m};
        idx = isnan(x) | isnan(y);
        x(idx) = [];
        y(idx) = [];

        all_x = [all_x; x(:)];
        all_y = [all_y; y(:)];
        all_c = [all_c; repmat(kk, length(x), 1)];  % 记录颜色组别

        % 绘制当前组的散点图
        scatter(x, y, 25, 'filled', 'MarkerFaceColor', colormap(kk,:), 'MarkerEdgeColor', 'none');
    end

    % === 进行合并后的拟合 ===
    if numel(all_x) >= 2  % 至少有两点才能拟合
        ft = fittype('poly1');  % 线性拟合
        [f, gof] = fit(all_x, all_y, ft);

        % 拟合曲线
        x_fit = linspace(min(all_x)-5, max(all_x)+5, 100);
        y_fit = feval(f, x_fit);
        plot(x_fit, y_fit, 'k-', 'LineWidth', 2);  % 黑色拟合线

        % 相关性分析
        [r, pValue] = corrcoef(all_x, all_y);
        rValue = r(1,2);
        pValue = pValue(1,2);
        text(min(all_x), max(all_y), sprintf('r = %.2f\\np = %.4f', rValue, pValue), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 8);
    end

    xlim([0, 120])
    ylabel(ParamName{m})
    set(gca, 'PositionConstraint', 'innerposition')
    box off
end



