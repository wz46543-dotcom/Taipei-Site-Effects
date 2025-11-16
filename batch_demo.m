% clear all;DataPath='F:\Job\TRF_KOsmth\Data\';OutputPath='F:\Job\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=10;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.,0);

%clear all;DataPath='F:\Job\TRF_KOsmth\Data\';OutputPath='F:\Job\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=0;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.75,1);

%clear all;DataPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\';OutputPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\';Ti=0.;Tf=40.;Tei=10;NP=4000;input_TRF;RnC_TRF('output','input','surface','bedrock','ground_response',0.,2);

clear all;DataPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Data\tainan_case\2016-02-05\';OutputPath='D:\signal_program\TRF_KOsmth\TRF_KOsmth\Analysis\tainan_case\2016-02-05\';
Ti=0.;Tf=120.;Tei=35;NP=4386;input_TRF;
F_min=0.2;F_max=30;RnC_TRF('N_2016-02-05CHY078.txt','U_2016-02-05CHY078.txt','N','U','2016-02-05CHY078',0.,2);

 % clear;DataPath='C:\Users\Asus\Desktop\input\Magnitude smaller than 5\2014-02-11\';OutputPath='C:\Users\Asus\Desktop\output\TRF\5以下\2014-02-11\';
 % Ti=0;Tf=120;Tei=66.5;NP=2000;input_TRF;
 % maker_TRF('E_TAP014','N_TAP014','U_TAP014','E','N','U','TAP014_2014-02-11',0,2.);