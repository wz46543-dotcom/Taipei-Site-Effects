function XQ_TRF(RE,RN,RU,nt,dt,LabE,LabN,LabU,ps,overlap,itp)
% XQ_TRF (Modified: Custom Excel Headers & Smoothed Data First)
%
% Updates:
% 1. Filename format: TR_combined_(E and N)wrt(U)_Station_Date.xlsx
% 2. Column Order: f, Smoothed(E), Smoothed(N), Raw(E), Raw(N)
% 3. Headers: Custom strings like 'TR-sm-E/U(f)'

global F_min F_max  T_i T_f T_ei T_ef NN
global T_tick F_tick 
global Data_Path Output_Path

T_unit='sec';
R_unit='gal'; 

% --- 檢查時間窗 ---
if (T_ei+NN*dt-dt > T_f)
    disp('Error: T_ef > T_f (Window exceeds file duration)')
    return
end
if (T_ei < T_i)
    disp('Error: T_ei < T_i')
    return
end

N = nt;
M = fix(N/2);

% --- 1. 全段訊號 FFT ---
RfE = fft(RE,N)*dt;
RfN = fft(RN,N)*dt;
RfU = fft(RU,N)*dt;

df = 1/(N*dt);
f = [0:df:M*df -(M-1)*df:df:-df];

% --- 2. 濾波 (Filter) ---
czero = complex(0,0);

% High pass
kf = find(abs(f)<F_min);
RfE(kf)=czero; RfN(kf)=czero; RfU(kf)=czero;

% Low pass
if (F_max > 1./2./dt | F_max == 0.)
    F_max = 1./2./dt;
end
kf = find(abs(f)>F_max);
RfE(kf)=czero; RfN(kf)=czero; RfU(kf)=czero;

% 反轉換回時間域
RE_fil = real(ifft(RfE,N)/dt); RE_fil = RE_fil(1:N);
RN_fil = real(ifft(RfN,N)/dt); RN_fil = RN_fil(1:N);
RU_fil = real(ifft(RfU,N)/dt); RU_fil = RU_fil(1:N);

t = 0:dt:(N-1)*dt;

% 頻率軸準備
f = 0:df:M*df;
jf = find(f<=F_max);

% --- 3. 迴圈計算 ---
icount = 1;
while(T_ei+NN*dt-dt <= T_f)

    idx_start = fix((T_ei)/dt)+1;
    idx_end = idx_start + NN - 1;
    
    ReffE = RE_fil(idx_start:idx_end);
    ReffN = RN_fil(idx_start:idx_end);
    ReffU = RU_fil(idx_start:idx_end);

    T_ef = T_ei + NN*dt;
    T_effect = T_ef - T_ei;

    % --- 加窗 (Hamming Window) ---
    N_wd = NN;
    w = hamming(N_wd); 
    
    ReffE_wd = ReffE .* w';
    ReffN_wd = ReffN .* w';
    ReffU_wd = ReffU .* w';

    % --- 補零 (Zero Padding) ---
    N_zp_tail = NN*0;  
    ReffE_zp = [ReffE_wd zeros(1,N_zp_tail)];
    ReffN_zp = [ReffN_wd zeros(1,N_zp_tail)];
    ReffU_zp = [ReffU_wd zeros(1,N_zp_tail)];

    N_zp = length(ReffE_zp);
    
    RfeE_zp = fft(ReffE_zp,N_zp)*dt;
    RfeN_zp = fft(ReffN_zp,N_zp)*dt;
    RfeU_zp = fft(ReffU_zp,N_zp)*dt;

    M_zp = fix(N_zp/2);
    df_zp = 1/(N_zp*dt);
    f_zp = 0:df_zp:M_zp*df_zp;

    % 截取有效頻率範圍
    jf_zp = find(f_zp<=F_max);
    f_zp = f_zp(jf_zp);
    RfeE_zp = RfeE_zp(jf_zp);
    RfeN_zp = RfeN_zp(jf_zp);
    RfeU_zp = RfeU_zp(jf_zp);

    % --- 4. 計算頻譜比 (Spectral Ratio) - 原始數據 ---
    TR_E = RfeE_zp ./ RfeU_zp; % E/U Raw
    TR_N = RfeN_zp ./ RfeU_zp; % N/U Raw

    czero = complex(0,0);
    kf_zp = find(f_zp < F_min);
    TR_E(kf_zp) = czero;
    TR_N(kf_zp) = czero;

    % --- 5. 平滑化 (Smoothing) - 平滑數據 ---
    if (itp == 2)
        TR_E_sm = smoothSpectra(TR_E,'w',20,'method','konno-ohmachi','b',20','debug','False');
        TR_N_sm = smoothSpectra(TR_N,'w',20,'method','konno-ohmachi','b',20','debug','False');
        
        % 抓取 E/U 峰值 (僅用於標示圖片)
        range_idx = fix(F_min/df_zp)+1:fix(10/df_zp); 
        if isempty(range_idx), range_idx = 1:length(f_zp); end
        
        [TRmax_E, maxfID_E] = max(abs(TR_E_sm(range_idx)));
        maxfID_E = maxfID_E + range_idx(1) - 1; 
        f_pre_E = f_zp(maxfID_E);
        TR_sm_pre_E = abs(TRmax_E);

        % 抓取 N/U 峰值 (僅用於標示圖片)
        [TRmax_N, maxfID_N] = max(abs(TR_N_sm(range_idx)));
        maxfID_N = maxfID_N + range_idx(1) - 1;
        f_pre_N = f_zp(maxfID_N);
        TR_sm_pre_N = abs(TRmax_N);
    end

    % ==========================================================
    % --- 6. 繪圖 Figure 1: E w.r.t U ---
    % ==========================================================
    figure(1); clf;
    set(gcf, 'Name', ['TR_(' LabE ')wrt(' LabU ')'], 'NumberTitle', 'off');
    
    subplot(3,1,1)
    plot(t, RE_fil, 'b-', 'LineWidth', 0.2);
    clear TITLE
    TITLE{1} = ['Remark: ' ps ];
    TITLE{2} = ['Time Window:[ ' num2str(T_ei,'%5.2f') ', ' num2str(T_ef,'%5.2f') '] sec ,  T_eff=' num2str(T_effect) ' sec,  Filter: [' num2str(F_min) ',' num2str(F_max) '] Hz'];
    title(TITLE, 'FontSize', 12, 'Interpreter', 'none', 'FontWeight', 'bold');
    ylabel(['a(t) (' R_unit ')'], 'FontSize', 10);
    set_time_axis(gca, T_tick, T_ei, T_ef);
    legend(['(' LabE ')'], 'Location', 'NorthEast');

    subplot(3,1,2)
    plot(t, RU_fil, 'b-', 'LineWidth', 0.2);
    ylabel(['a(t) (' R_unit ')'], 'FontSize', 10);
    xlabel(['Time (' T_unit ')'], 'FontSize', 10);
    set_time_axis(gca, T_tick, T_ei, T_ef);
    legend(['(' LabU ')'], 'Location', 'NorthEast');

    subplot(3,1,3)
    % 這裡先畫 Raw (紅色細線)
    semilogy(f_zp, abs(TR_E), 'r-', 'LineWidth', 0.3); hold on;
    if (itp == 2)
        % 再畫 Smoothed (綠色粗線)
        plot(f_zp, abs(TR_E_sm), 'g-', 'LineWidth', 2);
        text(f_pre_E, TR_sm_pre_E*1.1, ...
            ['(' num2str(f_pre_E,'%5.3f') ', ' num2str(TR_sm_pre_E,'%5.3f') ')\rightarrow '], ...
            'HorizontalAlignment', 'right', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 11);
    end
    hold off;
    
    xlim([0.2 10]); 
    ylim([0.01 100]); 
    
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Spectral Ratio', 'FontSize', 11);
    grid on;
    
    if (itp == 2)
        leg = legend(['(' LabE ') w.r.t (' LabU ')'], 'Smoothed');
    else
        leg = legend(['(' LabE ') w.r.t (' LabU ')']);
    end
    set(leg, 'FontSize', 9, 'Location', 'NorthEast');
    set(gcf, 'Position', [100, 50, 800, 700]);

    % 輸出圖片
    if (itp == 2)
        fig_out = [Output_Path 'TR_(' LabE ')wrt(' LabU ')_' ps '.png'];
        saveas(gcf, fig_out, 'png');
    end

    % ==========================================================
    % --- 7. 繪圖 Figure 2: N w.r.t U ---
    % ==========================================================
    figure(2); clf;
    set(gcf, 'Name', ['TR_(' LabN ')wrt(' LabU ')'], 'NumberTitle', 'off');
    
    subplot(3,1,1)
    plot(t, RN_fil, 'b-', 'LineWidth', 0.2);
    clear TITLE
    TITLE{1} = ['Remark: ' ps ];
    TITLE{2} = ['Time Window:[ ' num2str(T_ei,'%5.2f') ', ' num2str(T_ef,'%5.2f') '] sec ,  T_eff=' num2str(T_effect) ' sec,  Filter: [' num2str(F_min) ',' num2str(F_max) '] Hz'];
    title(TITLE, 'FontSize', 12, 'Interpreter', 'none', 'FontWeight', 'bold');
    ylabel(['a(t) (' R_unit ')'], 'FontSize', 10);
    set_time_axis(gca, T_tick, T_ei, T_ef);
    legend(['(' LabN ')'], 'Location', 'NorthEast');

    subplot(3,1,2)
    plot(t, RU_fil, 'b-', 'LineWidth', 0.2);
    ylabel(['a(t) (' R_unit ')'], 'FontSize', 10);
    xlabel(['Time (' T_unit ')'], 'FontSize', 10);
    set_time_axis(gca, T_tick, T_ei, T_ef);
    legend(['(' LabU ')'], 'Location', 'NorthEast');

    subplot(3,1,3)
    % 先畫 Raw
    semilogy(f_zp, abs(TR_N), 'r-', 'LineWidth', 0.3); hold on;
    if (itp == 2)
        % 再畫 Smoothed
        plot(f_zp, abs(TR_N_sm), 'g-', 'LineWidth', 2);
        text(f_pre_N, TR_sm_pre_N*1.1, ...
            ['(' num2str(f_pre_N,'%5.3f') ', ' num2str(TR_sm_pre_N,'%5.3f') ')\rightarrow '], ...
            'HorizontalAlignment', 'right', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 11);
    end
    hold off;
    
    xlim([0.2 10]);
    ylim([0.01 100]); 

    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Spectral Ratio', 'FontSize', 11);
    grid on;
    
    if (itp == 2)
        leg = legend(['(' LabN ') w.r.t (' LabU ')'], 'Smoothed');
    else
        leg = legend(['(' LabN ') w.r.t (' LabU ')']);
    end
    set(leg, 'FontSize', 9, 'Location', 'NorthEast');
    set(gcf, 'Position', [950, 50, 800, 700]);

    % 輸出圖片
    if (itp == 2)
        fig_out = [Output_Path 'TR_(' LabN ')wrt(' LabU ')_' ps '.png'];
        saveas(gcf, fig_out, 'png');
    end
    
    % ==========================================================
    % --- 8. 輸出 Excel 檔案 (新格式) ---
    % ==========================================================
    if (itp == 2)
        % 檔名格式：TR_combined_(E and N)wrt(U)_站名_日期.xlsx
        % (ps 變數應包含 "StationName_Date")
        Excel_Filename = [Output_Path 'TR_combined_(E and N)wrt(U)_' ps '.xlsx'];
        
        % 準備資料 (平滑在前面)
        v_f = f_zp(:);
        v_EU_sm  = abs(TR_E_sm(:)); % 平滑 E/U
        v_NU_sm  = abs(TR_N_sm(:)); % 平滑 N/U
        v_EU_raw = abs(TR_E(:));    % 原始 E/U
        v_NU_raw = abs(TR_N(:));    % 原始 N/U
        
        DataMat = [v_f, v_EU_sm, v_NU_sm, v_EU_raw, v_NU_raw];
        
        % 嘗試使用 writecell 來寫入帶有特殊符號的標題
        try
            % 自訂標題
            Header = {'f', 'TR-sm-E/U(f)', 'TR-sm-N/U(f)', 'TR-Raw-E/U(f)', 'TR-Raw-N/U(f)'};
            % 合併標題與數據
            OutputCell = [Header; num2cell(DataMat)];
            % 寫入
            writecell(OutputCell, Excel_Filename);
            fprintf('Excel output saved (writecell): %s\n', Excel_Filename);
            
        catch
            % 如果 MATLAB 版本較舊不支援 writecell，退回使用 writetable
            % 變數名稱必須符合 MATLAB 規則 (底線取代減號與括號)
            warning('writecell failed or not supported. Falling back to writetable with standard headers.');
            
            T = table(v_f, v_EU_sm, v_NU_sm, v_EU_raw, v_NU_raw, ...
                'VariableNames', {'f', 'TR_sm_E_U_f', 'TR_sm_N_U_f', 'TR_Raw_E_U_f', 'TR_Raw_N_U_f'});
            writetable(T, Excel_Filename);
            fprintf('Excel output saved (writetable): %s\n', Excel_Filename);
        end
    end

    % --- 9. 更新迴圈條件 ---
    if (itp == 1)
        icount = icount + 1;
        T_ei = T_ei + NN*(1.-overlap)*dt;
    else
        T_ei = T_f; % 結束迴圈
    end

end % End While

end

% --- 子程式: 設定時間軸樣式 ---
function set_time_axis(ax, T_tick, T_ei, T_ef)
    if( exist('T_tick','var') )
        set(ax,'XLim',[T_tick(1) T_tick(end)],'XTick',T_tick);
    end
    YLIMa = get(ax,'YLim');
    if( YLIMa(1) ~= -YLIMa(2) )
        MaxYLIMa = max(abs(YLIMa));
        YLIMa(1) = -MaxYLIMa;
        YLIMa(2) = +MaxYLIMa;
    end
    set(ax,'YLim',YLIMa);
    hold(ax, 'on');
    plot(ax, [T_ei T_ei], YLIMa, 'g--', 'LineWidth', 1.5);
    plot(ax, [T_ef T_ef], YLIMa, 'g--', 'LineWidth', 1.5);
    hold(ax, 'off');
    grid(ax, 'on');
end