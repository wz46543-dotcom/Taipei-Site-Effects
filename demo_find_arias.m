clear; clc;

% 1. --------------------------------
%    請設定您在 batch_demo.m 中使用的路徑和檔案
% --------------------------------
Data_Path = 'D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\tainan_case\2016-02-05\';
filename = 'E_TAP003.txt'; % 您要分析的檔案

full_path = [Data_Path, filename];

% 2. --------------------------------
%    設定您想查詢的能量百分比 (如同您的簡報 )
% --------------------------------
percentages_to_find = [5, 75, 95]; 

fprintf('===== 開始分析檔案: %s =====\n', filename);

% 3. --------------------------------
%    迴圈呼叫 find_arias_time 函數
%    (請確保 find_arias_time.m 與此腳本在同一個資料夾)
% --------------------------------
try
    results = [];
    for i = 1:length(percentages_to_find)
        p = percentages_to_find(i);
        time_point = find_arias_time(full_path, p);
        
        fprintf('  -> %g%% 能量對應時間點 (T_eff): %f 秒\n', p, time_point);
        results(i) = time_point;
    end
    
    % 4. (選用) 計算 5%-95% 有效歷時 (Duration) [cite: 163, 164]
    duration_5_95 = results(3) - results(1);
    fprintf('--------------------------------------------------\n');
    fprintf('  => 5%%-95%% 有效歷時 (Duration): %f 秒\n', duration_5_95);
    

catch ME
    fprintf('*** 分析時發生錯誤: %s ***\n', ME.message);
end