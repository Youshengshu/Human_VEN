function [spkwaveform,ttspan,...
    HW, AHP, Threshold, Amp,inflection_idx,initialSlope,max_slope,min_slope,...
    HWidx,AHPidx, Thresholdidx,min_slope2idx] = AP_Statistic(data,peaktime,tspan,Goodidx)

% INPUT:
% data:       voltage trace (vector)
% tspan:      time vector corresponding to data
% peaktime:   estimated spike peak times
% Goodidx:    indices of spikes to analyze

%% Preprocess

% Time resolution
dt = tspan(end) - tspan(end-1);

% Extract spike-aligned waveforms
spktime = peaktime(Goodidx);
pretime = 0.003;    % 3 ms before spike peak
postime = 0.004;    % 4 ms after spike peak

spkwaveform = [];

for i = 1:length(spktime)
    idx = round(spktime(i)/dt - pretime/dt : spktime(i)/dt + postime/dt);
    ttspan = linspace(-pretime, postime, numel(idx));
    spkwaveform = [spkwaveform, data(idx)];
end

%% Initialize
nspk = length(spktime);
inflection_idx = nan(1, nspk);
initialSlope = nan(1, nspk);
Thresholdidx = nan(1, nspk);
Threshold = nan(1, nspk);
Amp = nan(1, nspk);
AHPidx = nan(1, nspk);
AHP = nan(1, nspk);
HWidx = nan(nspk,2);
HW = nan(1, nspk);
max_slope = nan(1, nspk);
min_slope = nan(1, nspk);
max_slopeidx = nan(1, nspk);
min_slopeidx = nan(1, nspk);
min_slope2idx = nan(nspk,2);

%% Analyze spike features
for i = 1:length(spktime)
    data = spkwaveform(:,i);
    slope = [nan; diff(data)/(dt*1000)];          % first derivative in V/s
    slope2 = [nan; diff(slope)/(dt*1000)];        % second derivative in V/s^2
    slope2 = smooth(slope2, 5);

     % Find 2 major peaks in second derivative
    timewindow = [10:round(pretime/dt)];
    min_slope2idx(i,:) = peakfinder(slope2(timewindow),4,200,1,1) + 10 - 1;
    peaks = min_slope2idx(i,:);
    peaks = peaks(~isnan(peaks));  % 清除 NaN
    num_peaks = numel(peaks);

    if num_peaks >= 2
        % 情况1：两个明显 peak
        idx1 = min(peaks(1:2));
        idx2 = max(peaks(1:2));
        search_range = idx1:idx2;
        [~, min_idx] = min(slope2(search_range));
        inflection_idx(i) = search_range(min_idx);
    else
        % 情况2：只有一个或没有 peak，使用 slope2 上升相
        if ~isempty(peaks)
            main_peak_idx = peaks(1);
        else
            main_peak_idx = round(pretime/dt);  % 默认窗口
        end
        % 在 peak 前的 30% ~ 100% 范围内寻找局部最低点
        search_range = max(2, round(main_peak_idx * 0.3)) : main_peak_idx;
        [~, min_idx] = min(slope2(search_range));
        inflection_idx(i) = search_range(min_idx);
    end


    % === 对应的初始斜率 ===
    initialSlope(i) = slope(inflection_idx(i));

    %% Threshold and other AP features
    idxA = find(slope(:) > 20);
    idxB = find(data(:) > -60);
    idxC = intersect(idxA(:), idxB(:));

    if length(spktime) > 1 && i < length(peaktime)
        endtime = (peaktime(i+1) - spktime(i));
    else
        endtime = postime * 0.3;
    end

    if ~isempty(idxC)
        Thresholdidx(i) = idxC(1);
        Threshold(i) = data(Thresholdidx(i));
        AHPidx(i) = find(data(ttspan > 0 & ttspan < endtime) == min(data(ttspan > 0 & ttspan < endtime)), 1, 'first') + sum(ttspan < 0);
        AHP(i) = data(Thresholdidx(i)) - data(AHPidx(i));
        Amp(i) = max(data) - data(Thresholdidx(i));
        try
            HWidx(i,:) = eventedgefinder(data, data(Thresholdidx(i)) + Amp(i)/2, 1/dt, 0.0001, 0.0001, 1, 0);
        catch
            HWidx(i,:) = [nan nan];
        end
        HW(i) = (HWidx(i,2) - HWidx(i,1)) * dt * 1000;

        if ~isnan(HWidx(i,2))
            [~, AHPidx(i)] = min(data(HWidx(i,2):end-50));
            AHPidx(i) = AHPidx(i) + HWidx(i,2);
            AHP(i) = data(Thresholdidx(i)) - data(AHPidx(i));

            [~, max_slopeidx(i)] = max(slope(Thresholdidx(i):AHPidx(i)));
            max_slopeidx(i) = max_slopeidx(i) + Thresholdidx(i);
            max_slope(i) = slope(max_slopeidx(i));

            [~, min_slopeidx(i)] = min(slope(Thresholdidx(i):AHPidx(i)));
            min_slopeidx(i) = min_slopeidx(i) + Thresholdidx(i);
            min_slope(i) = slope(min_slopeidx(i));
        end
    end
end