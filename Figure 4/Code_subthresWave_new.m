clear all; clc
close all;

% 初始化数据结构用于保存统计结果
statsResults = struct();
statsCounter = 1;

Fidx{1} = dir('DATA_FIcurve*ETPC*.mat');
Fidx{2} = dir('DATA_FIcurve*cVEN*.mat');

colormap = [0.12,0.8,1; 1,0,0; 0.96,0.64,0.38; 0.75,0.75,0.75];  % 深蓝 红色 浅黄 浅灰
filename_Record = [];
subwave_100pA_Record = [];

% 数据收集和处理
for kk = 1:length(Fidx)
    for k = 1:length(Fidx{kk})
        tempVmFilter_Record = [];
        StimAmp = [];
        %% loading data
        filename = Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k), '/', num2str(length(Fidx{kk})), '--', num2str(kk), '/', num2str(length(Fidx)), '--', filename])
        dtV = tspanV(2) - tspanV(1);
        
        figure(1), clf
        for m = 1:size(waveformV, 2)
            tempVm = waveformV(:, m);
            windowSize = 1500;
            tempVmFiltered = medfilt1(tempVm, windowSize);
            tempVmFilter_Record = [tempVmFilter_Record, tempVmFiltered];
        end
        
        subwave_Record{kk}{k} = tempVmFilter_Record;     
        stimAmp_Record{kk}{k} = StimAmp;
    end
end

%% Bins 配置
pretime = 0.1;
postime = 0.8;

Binwidth = 20;
BinL = 0:Binwidth:180;
BinR = BinL + Binwidth;
BinsFIcurveS = [BinL', BinR'];

% 使用相同的Bin配置
Bins = BinsFIcurveS;

figure(2), clf
for kk = 1:length(Fidx)
    for k = 1:length(Fidx{kk})
        filename = Fidx{kk}(k).name;
        disp([num2str(k), '/', num2str(length(Fidx{kk})), '--', num2str(kk), '/', num2str(length(Fidx)), '--', filename])
        
        dataidx = round(pretime/dtV : (pretime+0.5)/dtV);
        baseidx = round(1 : (pretime-0.05)/dtV);
        Area = trapz(subwave_Record{kk}{k}(dataidx, :) - mean(subwave_Record{kk}{k}(baseidx, :)));
        [current, iidx] = sort(stimAmp_Record{kk}{k}, 'ascend');
        Area = Area(iidx);
        
        [current_unique, ~, ~] = unique(current);
        Area_unique = Area;
        
        current_color = colormap(kk, :);        
        color_with_alpha = [current_color, 0.4];
        
        plot(current_unique(current_unique < 200), Area_unique(current_unique < 200), 'Color', color_with_alpha);
        hold on;
        
        for i = 1:size(Bins, 1)
            idx = find(current_unique >= Bins(i, 1) & current_unique < Bins(i, 2));
            if ~isempty(idx)
                subArea_Bin{kk}(k, i) = mean(Area_unique(idx));
            else
                subArea_Bin{kk}(k, i) = nan;
            end
        end            
    end
end

%% 拟合和统计检验
figure(2), hold on
for kk = 1:length(Fidx)
    n = 1:length(Bins);
    x = mean(Bins, 2);
    data{kk} = subArea_Bin{kk}(:, n);
    data{kk}(find(data{kk} == 0)) = nan;
    y{kk} = nanmean(data{kk}, 1);
    N = sum(~isnan(data{kk}), 1);
    sem = nanstd(data{kk}, [], 1) ./ sqrt(N);
    
    idx = find(~isnan(y{kk}));
    x = x(idx);
    y{kk} = y{kk}(idx);
    data{kk} = data{kk}(:, idx);
    sem = sem(idx);
    xx = linspace(x(1), x(end), 100);
    
    [ModelLogn, gof] = fit(x(:), y{kk}(:), 'poly2');
    yy = feval(ModelLogn, xx);    
    plot(xx, yy, '-', 'color', colormap(kk, :));
    hold on
    errorbar(x, y{kk}, sem, 'ro', 'color', colormap(kk, :), 'markerfacecolor', colormap(kk, :), 'markersize', 8)
end
ylabel('Area')
xlabel('Current injection (pA)')
box off

%% ======================= 修改后的统计检验部分 =======================
% 创建表格存储结果
resultsTable = table();
binCenters = mean(Bins, 2);

% 确保有swtest函数（正态性检验）
if ~exist('swtest', 'file')
    error('缺少swtest函数，请确保其在MATLAB路径中');
end

% 为每个Bin添加统计检验
for binIdx = 1:size(Bins, 1)
    group1_data = subArea_Bin{1}(:, binIdx);
    group2_data = subArea_Bin{2}(:, binIdx);
    
    % 移除NaN值
    group1_data = group1_data(~isnan(group1_data));
    group2_data = group2_data(~isnan(group2_data));
    
    % 检查是否有足够的数据
    if numel(group1_data) < 4 || numel(group2_data) < 4
        fprintf('Bin %d (%d-%d pA): 数据不足，跳过统计检验\n', binIdx, Bins(binIdx, 1), Bins(binIdx, 2));
        testMethod = '数据不足';
        pValue = NaN;
        continue;
    end
    
    % 初始化变量
    testMethod = '';
    pValue = NaN;
    
    % 按照您提供的逻辑进行检验
    if ~isempty(group1_data) && ~isempty(group2_data)
        % Shapiro-Wilk正态性检验 (α=0.05)
        normality1 = swtest(group1_data, 0.05);
        normality2 = swtest(group2_data, 0.05);
        
        % 两组都符合正态分布 → 双样本t检验
        if normality1 == 0 && normality2 == 0
            [~, pValue] = ttest2(group1_data, group2_data);
            testMethod = 'Two samples t-test';
            
        % 至少一组不符合正态分布 → Mann-Whitney U检验
        else
            [pValue, ~] = ranksum(group1_data, group2_data);
            testMethod = 'Wilcoxon rank sum test';
        end
        
        % 计算效应量 Cohen's d (无论使用何种检验都计算)
        n1 = length(group1_data);
        n2 = length(group2_data);
        pooled_std = sqrt(((n1-1)*var(group1_data) + (n2-1)*var(group2_data)) / (n1 + n2 - 2));
        cohen_d = (mean(group1_data) - mean(group2_data)) / pooled_std;
    end
    
    % 保存结果
    resultsTable.BinStart(binIdx) = Bins(binIdx, 1);
    resultsTable.BinEnd(binIdx) = Bins(binIdx, 2);
    resultsTable.n_Group1(binIdx) = n1;
    resultsTable.n_Group2(binIdx) = n2;
    resultsTable.mean_Group1(binIdx) = mean(group1_data);
    resultsTable.mean_Group2(binIdx) = mean(group2_data);
    resultsTable.std_Group1(binIdx) = std(group1_data);
    resultsTable.std_Group2(binIdx) = std(group2_data);
    resultsTable.TestMethod(binIdx) = {testMethod};
    resultsTable.pValue(binIdx) = pValue;
    resultsTable.Cohen_d(binIdx) = cohen_d;
    
    % 打印结果
    fprintf('\nBin %d (%d-%d pA) 统计结果:\n', binIdx, Bins(binIdx, 1), Bins(binIdx, 2));
    fprintf('-----------------------------------\n');
    fprintf('检验方法: %s\n', testMethod);
    fprintf('p值 = %.4f\n', pValue);
    fprintf('Cohen''s d = %.2f\n', cohen_d);
    fprintf('Group1: n=%d, mean=%.2f±%.2f\n', n1, mean(group1_data), std(group1_data));
    fprintf('Group2: n=%d, mean=%.2f±%.2f\n', n2, mean(group2_data), std(group2_data));
end

% 多重比较校正 (FDR)
pValues = resultsTable.pValue;
[~, ~, ~, adj_p] = fdr_bh(pValues, 0.05, 'pdep');

% 添加校正后的p值到表格
resultsTable.Adj_pValue = adj_p;

% 保存结果到CSV文件
writetable(resultsTable, 'StatisticalResults.csv');
save('StatisticalResults.mat', 'resultsTable');

fprintf('\n统计结果已保存到 StatisticalResults.csv 和 StatisticalResults.mat\n');

%% 辅助函数: FDR校正 (与之前相同)
function [h, crit_p, adj_p, index] = fdr_bh(pvals, q, method)
% FDR_BH Benjamini-Hochberg false discovery rate procedure
    if nargin < 3
        method = 'pdep';
    end

    pvals = pvals(:);
    m = length(pvals);
    
    % Sort p-values and remember original indices
    [sorted_p, index] = sort(pvals);
    
    % Adjust p-values
    if strcmpi(method, 'pdep')
        adj_p = m * sorted_p ./ (1:m)';
    elseif strcmpi(method, 'dep')
        c = sum(1./(1:m));
        adj_p = m * sorted_p * c ./ (1:m)';
    else
        error('Unknown method');
    end
    
    % Ensure adjusted p-values don't decrease
    adj_p = cummin(adj_p, 'reverse');
    
    % Find critical p-value
    crit_idx = find(sorted_p <= (1:m)' * q / m, 1, 'last');
    if isempty(crit_idx)
        crit_p = 0;
    else
        crit_p = sorted_p(crit_idx);
    end
    
    % Determine rejected hypotheses
    h = pvals <= crit_p;
    
    % Restore original order for adj_p
    [~, orig_idx] = sort(index);
    adj_p = adj_p(orig_idx);
end