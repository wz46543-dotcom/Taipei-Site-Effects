% clear all;DataPath='F:\Job\TRF_KOsmth\Data\';OutputPath='F:\Job\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=10;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.,0);

%clear all;DataPath='F:\Job\TRF_KOsmth\Data\';OutputPath='F:\Job\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=0;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.75,1);

%clear all;DataPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\';OutputPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=10;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.,2);

% clear all;DataPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\tainan_case\2016-02-05\';OutputPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\tainan_case\2016-02-05\';
% Ti=0.;Tf=120.;Tei=35;NP=4386;input_TRF;
% F_min=0.2;F_max=30;RnC_TRF('E_2016-02-05CHY078.txt','U_2016-02-05CHY078.txt','E','U','2016-02-05CHY078',0.,2);

 clear all;
% 設定資料路徑與輸出路徑 (請依您的電腦修改)
DataPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\tainan_case\2010-03-04\';
OutputPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\tainan_case\2010-03-04\';


% 設定全域變數
Ti=0.; Tf=120.; Tei=10.00; NP=4680;
input_TRF;

% 頻率範圍設定
F_min=0.2; F_max=30;

% 修改後的呼叫方式: 傳入 (E檔名, N檔名, U檔名, 'E', 'N', 'U', 備註, overlap, 模式)
% 模式 2 代表進行 Konno-Ohmachi 平滑化
RnC_TRF('E_2010-03-04CHY021.txt','N_2010-03-04CHY021.txt','U_2010-03-04CHY021.txt','E', 'N', 'U', '2010-03-04CHY021', 0., 2);
