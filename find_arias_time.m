function time_point = find_arias_time(full_filepath, percentage)
% =========================================================================
% find_arias_time: 根據愛氏震度 (Arias Intensity) 計算指定能量百分比的時間點
%
% [輸入]
%   full_filepath: 完整檔案路徑 (包含檔名)，例如 'D:\...Data\N_...CHY078.txt'
%   percentage:    要查詢的能量百分比 (請輸入數字，例如 95)
%
% [輸出]
%   time_point:    計算出的時間點 (秒)
% =========================================================================

% --- 1. 讀取地震訊號檔案 ---
% 這個讀取邏輯與您提供的 RnC_TRF.m 完全相同
[fid, message] = fopen(full_filepath, 'r');
if fid == -1
    error('檔案開啟失敗: %s\n%s', full_filepath, message);
end

mult = fscanf(fid, '%g', 1); % 讀取第1行 (Multiplier)
dt = fscanf(fid, '%g', 1);   % 讀取第2行 (Time step, e.g., 0.005)
[accel_data, count] = fscanf(fid, '%g', [1, inf]); % 讀取第3行之後所有資料
accel_data = accel_data * mult; % 套用乘積
fclose(fid);

if isempty(accel_data)
    error('檔案中沒有讀取到地震訊號資料: %s', full_filepath);
end

% --- 2. 計算愛氏震度 (Arias Intensity) 的累計能量 ---
% 愛氏震度 Ia = (pi/2g) * integral(a(t)^2) dt
% 我們只需要計算「歸一化」的能量，所以 (pi/2g) 常數可以省略

% (a) 計算 a(t)^2
a_squared = accel_data .^ 2;

% (b) 計算累計總和 (即 integral(a^2) 的離散版本)
%     使用 cumsum (Cumulative Sum)
E_cumulative = cumsum(a_squared);

% (c) 取得總能量 (最後一個值)
E_total = E_cumulative(end);

% (d) 歸一化，得到 0 到 1 之間的累計能量百分比曲線
E_normalized = E_cumulative / E_total;

% --- 3. 找出目標百分比對應的時間點 ---
target_ratio = percentage / 100.0; % 將 95 轉換為 0.95

% (a) 尋找第一個 E_normalized >= target_ratio 的「索引」(index)
index = find(E_normalized >= target_ratio, 1, 'first');

if isempty(index)
    % 如果找不到 (例如輸入 110%)，則回傳總時間
    warning('無法找到 %g%% 的能量點。回傳總歷時。', percentage);
    time_point = (count - 1) * dt;
    return;
end

% (b) 根據索引和時間間隔 (dt) 計算時間
% (索引 1 對應 t=0, 索引 2 對應 t=dt, ... 索引 i 對應 t=(i-1)*dt)
time_point = (index - 1) * dt;

end