%% 地震資料自動化分析 (變數名稱修正版)
% 1. 修正全域變數名稱 (Ti->T_i, Tf->T_f, Tei->T_ei) 以匹配 XQ_TRF.m
% 2. 包含之前所有的修正 (自動對齊長度、檔名簡化、資料夾美化)
%#ok<*GVMIS> 
%#ok<*NUSED>
%#ok<*UNRCH>

clear; clc; close all;
warning('off','all');

% ==============================================================================
% 1. 路徑設定
% ==============================================================================
Raw_Root = 'D:\seismo\out_20-25TXT\ML1-'; 
Data_Root = 'D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\ML6-';
Anal_Root = 'D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\ML6-';

Start_P = 5; End_P = 95;    
F_min_v = 0.2; F_max_v = 30;

% ==============================================================================
% 2. 執行確認
% ==============================================================================
if ~exist(Raw_Root, 'dir')
    Raw_Root = uigetdir('D:\', '請選擇包含地震資料的根目錄');
    if Raw_Root == 0, error('未選擇資料夾，程式結束'); end
end

E_list = dir(Raw_Root);
E_list = E_list([E_list.isdir]); 

fprintf('\n=== 開始執行 (全域變數名稱已修正: T_f, T_ei...) ===\n來源路徑: %s\n', Raw_Root);

for i = 1:length(E_list)
    if strcmp(E_list(i).name, '.') || strcmp(E_list(i).name, '..'), continue; end
    
    Ev_Name_Source = E_list(i).name; 
    Ev_Path = fullfile(Raw_Root, Ev_Name_Source); 
    
    % 資料夾名稱美化 (YYYY_MM_DD_HH_MM_SS)
    if length(Ev_Name_Source) == 14
        Ev_Name_Out = sprintf('%s_%s_%s_%s_%s_%s', ...
            Ev_Name_Source(1:4), Ev_Name_Source(5:6), Ev_Name_Source(7:8), ...
            Ev_Name_Source(9:10), Ev_Name_Source(11:12), Ev_Name_Source(13:14));
    else
        Ev_Name_Out = Ev_Name_Source; 
    end
    
    S_list = dir(fullfile(Ev_Path, 'TW.*'));
    S_list = S_list([S_list.isdir]);
    
    if isempty(S_list), continue; end
    
    fprintf('\n>>> 進入地震: %s -> 輸出為: %s\n', Ev_Name_Source, Ev_Name_Out);
    
    % 建立輸出目錄
    Ev_Data_P = fullfile(Data_Root, Ev_Name_Out); 
    Ev_Anal_P = fullfile(Anal_Root, Ev_Name_Out);
    if ~exist(Ev_Data_P, 'dir'), mkdir(Ev_Data_P); end
    if ~exist(Ev_Anal_P, 'dir'), mkdir(Ev_Anal_P); end

    for j = 1:length(S_list)
        S_Folder = S_list(j).name;
        S_ID = strrep(S_Folder, 'TW.', ''); 
        Full_S_Path = fullfile(Ev_Path, S_Folder);
        
        % --- A. 讀取與對齊 ---
        data_E = []; dt_E = 0;
        data_N = []; dt_N = 0;
        data_U = []; dt_U = 0;
        
        [data_E, dt_E] = read_component_data(Full_S_Path, {'HLE'});
        [data_N, dt_N] = read_component_data(Full_S_Path, {'HLN'});
        [data_U, dt_U] = read_component_data(Full_S_Path, {'HLU', 'HLZ'}); 
        
        if isempty(data_E) || isempty(data_N) || isempty(data_U)
            fprintf('    [跳過] 測站 %s 資料不全\n', S_ID);
            continue;
        end
        
        % 檢查取樣率
        if abs(dt_E - dt_N) > 1e-6
             fprintf('    [跳過] 測站 %s 取樣率不一致\n', S_ID);
             continue;
        end
        dt = dt_E;
        
        % 自動對齊長度
        len_E = length(data_E); len_N = length(data_N); len_U = length(data_U);
        min_len = min([len_E, len_N, len_U]);
        
        if (len_E ~= min_len || len_N ~= min_len || len_U ~= min_len)
            data_E = data_E(1:min_len);
            data_N = data_N(1:min_len);
            data_U = data_U(1:min_len);
        end
        
        % --- B. 寫入文字檔 (格式: 1, dt, data 分行) ---
        write_txt_file(Ev_Data_P, 'E', S_ID, dt, data_E);
        write_txt_file(Ev_Data_P, 'N', S_ID, dt, data_N);
        write_txt_file(Ev_Data_P, 'U', S_ID, dt, data_U);

        % --- C. 執行 HVSR 分析 ---
        try
            % === 關鍵修正：變數名稱必須與 XQ_TRF.m 完全一致 (加底線) ===
            global Data_Path Output_Path T_i T_f T_ei T_ef NN F_min F_max
            global T_tick F_tick
            
            Data_Path = [Ev_Data_P filesep];
            Output_Path = [Ev_Anal_P filesep];
            
            % 計算時間窗
            e_txt_full = fullfile(Ev_Data_P, sprintf('E_%s.txt', S_ID));
            
            % 計算起始與結束點
            t_start = find_arias_time(e_txt_full, Start_P); 
            t_end = find_arias_time(e_txt_full, End_P);
            
            % 設定到全域變數 (注意名稱是 T_ei 而不是 Tei)
            T_ei = t_start;
            NN = round((t_end - t_start) / dt);
            
            % 設定檔案總長度
            Total_Duration = min_len * dt; 
            T_f = Total_Duration; % 注意名稱是 T_f
            
            % 安全檢查：確保分析結束點不超過檔案長度
            if (T_ei + NN * dt) > Total_Duration
                NN = floor((Total_Duration - T_ei) / dt) - 5;
            end
            
            T_i = 0.; 
            F_min = F_min_v; F_max = F_max_v;
            
            % 設定繪圖刻度
            T_tick = 0:20:ceil(T_f); 
            F_tick = 0:2:30;

            Rem = sprintf('%s_%s', Ev_Name_Out, S_ID); 
            fnE = sprintf('E_%s.txt', S_ID);
            fnN = sprintf('N_%s.txt', S_ID);
            fnU = sprintf('U_%s.txt', S_ID);
            
            RnC_TRF(fnE, fnN, fnU, 'E', 'N', 'U', Rem, 0, 2);
            
            % --- D. 驗證產出 ---
            expected_excel = fullfile(Ev_Anal_P, sprintf('TR_combined_(E and N)wrt(U)_%s.xlsx', Rem));
            if exist(expected_excel, 'file')
                fprintf('    [成功] 測站 %s 分析完成\n', S_ID);
            else
                fprintf('    [失敗] 測站 %s 無產出。參數診斷: T_f=%.1f, T_ei=%.1f, NP=%d\n', ...
                    S_ID, T_f, T_ei, NN);
            end
            
        catch ME
            fprintf('    [崩潰] %s: %s\n', S_ID, ME.message);
        end
    end
end
fprintf('\n=== 全部作業結束，請檢查 Analysis 資料夾 ===\n');

% ==============================================================================
% 輔助函數
% ==============================================================================
function [data, dt] = read_component_data(folder_path, tags)
    data = []; dt = 0;
    all_files = dir(folder_path);
    target_file = '';
    
    for m = 1:length(all_files)
        fname = all_files(m).name;
        if all_files(m).isdir, continue; end
        if contains(fname, '.mseed', 'IgnoreCase', true), continue; end
        if contains(fname, '.png', 'IgnoreCase', true), continue; end
        
        for t = 1:length(tags)
            if contains(fname, tags{t}, 'IgnoreCase', true)
                target_file = fullfile(folder_path, fname);
                break;
            end
        end
        if ~isempty(target_file), break; end
    end
    
    if isempty(target_file), return; end
    
    fid = fopen(target_file, 'r');
    lines_read = 0;
    while ~feof(fid) && lines_read < 20
        tline = fgetl(fid);
        lines_read = lines_read + 1;
        if contains(tline, 'sps')
            tk = regexp(tline, '(\d+)\s+sps', 'tokens');
            if ~isempty(tk)
                dt = 1 / str2double(tk{1}{1});
                dc = textscan(fid, '%s %f');
                data = dc{2};
                if isrow(data), data = data'; end
                break;
            end
        end
    end
    fclose(fid);
end

function write_txt_file(out_path, comp, station, dt, data)
    out_f = sprintf('%s_%s.txt', comp, station);
    fid_o = fopen(fullfile(out_path, out_f), 'w');
    fprintf(fid_o, '1\n');           
    fprintf(fid_o, '%g\n', dt);      
    fprintf(fid_o, '%.10e\n', data); 
    fclose(fid_o);
end