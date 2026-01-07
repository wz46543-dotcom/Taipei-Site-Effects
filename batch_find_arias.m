clear; clc;

% ==============================================================================
% 1. 使用者設定區
% ==============================================================================
Data_Path = 'D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\tainan_case\2010-03-04\';
StationName = '2010-03-04CHY021';             % 測站名稱
Tei_input = 10;                   % 您手動決定的 S 波起始時間 (Tei)
percentages_to_find = [5, 75, 95];  % 您想查詢的能量百分比 
% ==============================================================================

filename = ['E_' StationName '.txt']; 
full_path = [Data_Path, filename];
fprintf('===== 開始分析檔案: %s =====\n', filename);

try
    % --- 2. 讀取檔案取得 dt (為了計算 NP) ---
    [fid, message] = fopen(full_path, 'r');
    if fid == -1
        error('無法開啟檔案: %s \n錯誤訊息: %s', full_path, message);
    end
    fscanf(fid, '%g', 1);      % 跳過 mult
    dt = fscanf(fid, '%g', 1); % 讀取 dt
    fclose(fid);

    % --- 3. [功能 1] 顯示各百分比的時間點 (回復原功能) ---
    fprintf('\n--- [參考資訊] 愛氏震度能量時間點 ---\n');
    T_results = [];
    for i = 1:length(percentages_to_find)
        p = percentages_to_find(i);
        t_val = find_arias_time(full_path, p);
        fprintf('  -> %g%% 能量時間點: %.4f 秒\n', p, t_val);
        T_results(i) = t_val;
    end
    
    % 計算 5%-95% 有效歷時 (如果列表中有這兩個值)
    idx_5 = find(percentages_to_find == 5);
    idx_95 = find(percentages_to_find == 95);
    if ~isempty(idx_5) && ~isempty(idx_95)
        dur = T_results(idx_95) - T_results(idx_5);
        fprintf('  => 5%%-95%% 有效歷時 (Ds): %.4f 秒\n', dur);
    end

    % --- 4. [功能 2] 生成 Batch 程式碼 ---
    % 取得 95% 時間點 (作為結束點 T_ef)
    if ~isempty(idx_95)
        T_end = T_results(idx_95);
    else
        T_end = find_arias_time(full_path, 95); % 如果上面沒設95%，這裡補算
    end
    
    fprintf('\n--- [生成代碼] 參數計算 ---\n');
    fprintf('  -> 手動設定起始 (Tei): %.2f 秒\n', Tei_input);
    fprintf('  -> 自動計算結束 (95%%): %.4f 秒\n', T_end);

    if T_end <= Tei_input
        warning('注意: 95%% 結束時間 (%.2f) 小於 起始時間 (%.2f)，NP 計算可能錯誤。', T_end, Tei_input);
    end
    
    Duration_calc = T_end - Tei_input;
    NP = round(Duration_calc / dt); % 四捨五入取整數
    
    fprintf('  -> 計算結果: 擷取長度 = %.4f 秒, dt = %g, 點數 NP = %d\n', Duration_calc, dt, NP);

    % --- 5. 印出 Code Block ---
    fprintf('\n================ 請複製下方文字到 batch_demo.m ================\n\n');
    
    fprintf('%% 設定全域變數\n');
    fprintf('Ti=0.; Tf=120.; Tei=%.2f; NP=%d;\n', Tei_input, NP);
    fprintf('input_TRF;\n\n');
    
    fprintf('%% 頻率範圍設定\n');
    fprintf('F_min=0.2; F_max=30;\n\n');
    
    fprintf('%% 修改後的呼叫方式: 傳入 (E檔名, N檔名, U檔名, ''E'', ''N'', ''U'', 備註, overlap, 模式)\n');
    fprintf('%% 模式 2 代表進行 Konno-Ohmachi 平滑化\n');
    
    % 自動組合檔名字串
    str_E = sprintf('''E_%s.txt''', StationName);
    str_N = sprintf('''N_%s.txt''', StationName);
    str_U = sprintf('''U_%s.txt''', StationName);
    str_Rem = sprintf('''%s''', StationName); % 備註
    
    fprintf('RnC_TRF(%s,%s,%s,''E'', ''N'', ''U'', %s, 0., 2);\n', str_E, str_N, str_U, str_Rem);
    
    fprintf('\n==============================================================\n');

catch ME
    fprintf('*** 執行錯誤: %s ***\n', ME.message);
end