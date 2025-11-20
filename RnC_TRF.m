function RnC_TRF(fnameE, fnameN, fnameU, LabE, LabN, LabU, ps, overlap, itp)
% RnC_TRF - Read and Check Transfer Function (Modified for 3 Components)
% 修改版：同時讀取 E, N, U 三個分量的檔案，檢查格式後傳送給 XQ_TRF 進行計算。

% 宣告全域變數
global F_min F_max T_i T_f T_ei T_ef NN
global T_tick F_tick
global Data_Path Output_Path

% --- 1. 讀取 E 分量檔案 ---
F_nameE = [Data_Path fnameE];
[fidE, message] = fopen(F_nameE, 'r');
if fidE == -1
    disp(['Error opening E file: ' message]);
    return;
end
multE = fscanf(fidE, '%g', 1);      % 讀取倍率因子
dtE = fscanf(fidE, '%g', 1);        % 讀取取樣間隔
[RE, countE] = fscanf(fidE, '%g', [1, inf]); % 讀取數據
RE = RE * multE;
fclose(fidE);
ntE = countE;

% --- 2. 讀取 N 分量檔案 ---
F_nameN = [Data_Path fnameN];
[fidN, message] = fopen(F_nameN, 'r');
if fidN == -1
    disp(['Error opening N file: ' message]);
    return;
end
multN = fscanf(fidN, '%g', 1);
dtN = fscanf(fidN, '%g', 1);
[RN, countN] = fscanf(fidN, '%g', [1, inf]);
RN = RN * multN;
fclose(fidN);
ntN = countN;

% --- 3. 讀取 U 分量檔案 ---
F_nameU = [Data_Path fnameU];
[fidU, message] = fopen(F_nameU, 'r');
if fidU == -1
    disp(['Error opening U file: ' message]);
    return;
end
multU = fscanf(fidU, '%g', 1);
dtU = fscanf(fidU, '%g', 1);
[RU, countU] = fscanf(fidU, '%g', [1, inf]);
RU = RU * multU;
fclose(fidU);
ntU = countU;

% --- 4. 檢查相容性 (長度與取樣率) ---
if (ntE ~= ntU) || (ntN ~= ntU)
    disp('Error: Signal lengths do not match across E, N, U components.');
    return;
end

if (dtE ~= dtU) || (dtN ~= dtU)
    disp('Error: Sampling rates (dt) do not match across E, N, U components.');
    return;
end

% --- 5. 調整 T_f (若設定的結束時間超過檔案長度) ---
if (T_f > ntE * dtE)
    T_f = ntE * dtE;
end  

% --- 6. 呼叫核心計算程式 ---
XQ_TRF(RE, RN, RU, ntE, dtE, LabE, LabN, LabU, ps, overlap, itp);

% 這裡一定要有 end 來結束 function
end