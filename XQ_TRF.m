function XQ_TRF(R1,R2,nt,dt,CH1,CH2,ps,overlap,itp)
% % fname1: input file name of output signal 
% % fname2: input file name of input signal
% %         *output and input signals must have same length and dt. 
% % CH1: specified name of output signal 
% % CH2: specified name of input signal
% % ps: remarks
% % overlap (only effective for itp =1): overlap ratio of adjoining segments.
% % itp: index of execution mode. 
% %      if itp = 0, only calculate transfer function once [T_ei~(T_ei+NP*dt)]; 
% %      if itp = 1, calculate transfer functions of successive segments with a duration of NP*dt
% %                  starting from given T_ei, and T_ei shifts (1-overlap)*NP*dt for the next segment, until Tei+NP*dt > Tf 
% %      if itp = 2, same as itp = 0 but TRF is smoothed.

global F_min F_max  T_i T_f T_ei T_ef NN
global T_tick F_tick 
global Data_Path Output_Path

T_unit='sec';
R_unit='m/s^2';

if (T_ei+NN*dt-dt>T_f)
         disp('T_ef>T_f')
     return
end

if (T_ei < T_i)
         disp('T_ei<T_f')
     return
end

N=nt;
M=fix(N/2);
%
Rf1=fft(R1,N)*dt;
Rf2=fft(R2,N)*dt;
%disp(['R no.= ' num2str(length(R))]);
df=1/(N*dt);
f=[0:df:M*df -(M-1)*df:df:-df];

%
% HIGH pass filter
%
kf=find(abs(f)<F_min);
czero=complex(0,0);
Rf1(kf)=czero;
Rf2(kf)=czero;
%
% LOW pass filter
%
if (F_max > 1./2./dt | F_max == 0.)
    F_max = 1./2./dt;
end

kf=find(abs(f)>F_max);
Rf1(kf)=czero;
Rf2(kf)=czero;

R1_fil=real(ifft(Rf1,N)/dt);
R1_fil=R1_fil(1:N);
R2_fil=real(ifft(Rf2,N)/dt);
R2_fil=R2_fil(1:N);
t=0:dt:(N-1)*dt;

f=0:df:M*df;
jf=find(f<=F_max);
Rf1=Rf1(jf);
Rf2=Rf2(jf);
f=f(jf);

% if (itp==1)
%   Filename=[Output_Path ps  '_(' CH1 ')wrt(' CH2 ').pre'];
%   fid_pre=fopen(Filename,'w');
%   fprintf(fid_pre,'    t_cen     f_pre   TR_pre\n'); 
% end

%%+++++++++++++++++++++++++++++++++++++
icount=1;
while(T_ei+NN*dt-dt<=T_f)
%%+++++++++++++++++++++++++++++++++++++
Reff1=R1_fil(fix((T_ei)/dt)+1:fix((T_ei)/dt)+NN);
Reff2=R2_fil(fix((T_ei)/dt)+1:fix((T_ei)/dt)+NN);

T_ef=T_ei+NN*dt;
T_effect=T_ef-T_ei;
%disp([ 'Reff = ' , num2str(length(Reff))]);
%+++++++++++++++++++++++++++++++++++++
%Rfe1=fft(Reff1,NN)*dt;
%Rfe2=fft(Reff2,NN)*dt;
%clear Reff
%Rfe=Rfe(jf);

%%%
%%% ____________ windowing _________________
%%%
N_wd=NN;
M_wd=fix(N_wd/2);

j_wd=length(Reff1)-N_wd+1;
if j_wd<0
    j_wd=1;
end
t_wd=0:dt:(length(Reff1)-j_wd)*dt;
% w = hann(N_wd);
w = hamming(N_wd);
% w = ones(N_wd,1);
Reff1_wd=Reff1(j_wd:end);
Rh1_wd=Reff1_wd.*w';
Reff2_wd=Reff2(j_wd:end);
Rh2_wd=Reff2_wd.*w';

%% ____________ zero-padding _________________
N_zp_tail=NN*0;  %%%% length of zero-padding

Reff1_zp=[Rh1_wd zeros(1,N_zp_tail)];
Reff2_zp=[Rh2_wd zeros(1,N_zp_tail)];

N_zp=length(Reff1_zp);

Rfe1_zp=fft(Reff1_zp,N_zp)*dt;
Rfe2_zp=fft(Reff2_zp,N_zp)*dt;

M_zp=fix(N_zp/2);
df_zp=1/(N_zp*dt);
f_zp=0:df_zp:M_zp*df_zp;

jf_zp=find(f_zp<=F_max);
f_zp=f_zp(jf_zp);

Rfe1_zp=Rfe1_zp(jf_zp);

Rfe2_zp=Rfe2_zp(jf_zp);

% % % 
% % % =======  transfer function calculation  ===========
% % % 

TR=Rfe1_zp./Rfe2_zp;
TR=TR(jf_zp);

kf_zp=find(f_zp<F_min);
TR(kf_zp)=czero;

% % % =======  transfer function smoothing (itp==2)  ===========

if (itp==2)
    % TR_sm=smooth(TR,0.007,'rloess');
    % TR_sm=smoothdata(TR,'rloess',10);
    %TR_sm=kohmachi(TR,f_zp,20);
   
    TR_sm=smoothSpectra(TR,'w',20,'method','konno-ohmachi','b',20','debug','False');

end

if (icount==1)
  TR_ave=zeros(size(TR));
end 

TR_ave=TR_ave(jf_zp);

TR_ave=TR_ave+abs(TR);

if (itp==2)
  [TRmax,maxfID]=max(abs(TR_sm(fix(F_min/df_zp)+1:fix(10/df_zp))));
  TR_sm_pre=abs(TRmax);
  f_pre=f_zp(maxfID);

end

subplot(4,1,1)
plot(t,R1_fil,'b-','LineWidth',0.2);
clear TITLE
TITLE{1}=['Remark: ' ps ];
TITLE{2}=['Time Window:[ ' num2str(T_ei,'%5.2f') ', ' num2str(T_ef,'%5.2f') '] sec ,  T_eff=' num2str(T_effect) ' sec,  Filter: [' num2str(F_min) ',' num2str(F_max) '] Hz'];
title(TITLE,'FontSize',12,'Interpreter','none');
xlabel(['Time (' T_unit ')'],'FontSize',12);
clear YLABEL
YLABEL=['a(t) (' R_unit ')'];
ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
if( exist('T_tick','var') )
    set(gca,'XLim',[T_tick(1) T_tick(end)],'XTick',T_tick);
end
%axis([T_tick(1) T_tick(end) -0.005 0.005])
YLIMa=get(gca,'YLim');
if( YLIMa(1) ~= -YLIMa(2) )
  MaxYLIMa=max(abs(YLIMa));
  YLIMa(1)=-MaxYLIMa;
  YLIMa(2)=+MaxYLIMa;
end
set(gca,'YLim',YLIMa);
hold on
plot([T_ei T_ei],YLIMa,'g--','LineWidth',1.5);
plot([T_ef T_ef],YLIMa,'g--','LineWidth',1.5);
hold off
YTICKa=get(gca,'YTick');
if( length(YTICKa)==3 )
    YTICKa=YLIMa(1):0.25*(YLIMa(2)-YLIMa(1)):YLIMa(2);
end
set(gca,'YTick',YTICKa);
grid on
leg=legend(['(' CH1 ')']);
set(leg,'FontSize',12);

subplot(4,1,2)
plot(t,R2_fil,'b-','LineWidth',0.2);

xlabel(['Time (' T_unit ')'],'FontSize',12);
clear YLABEL
YLABEL=['a(t) (' R_unit ')'];
ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
if( exist('T_tick','var') )
    set(gca,'XLim',[T_tick(1) T_tick(end)],'XTick',T_tick);
end
%axis([T_tick(1) T_tick(end) -0.005 0.005])
YLIMa=get(gca,'YLim');
if( YLIMa(1) ~= -YLIMa(2) )
  MaxYLIMa=max(abs(YLIMa));
  YLIMa(1)=-MaxYLIMa;
  YLIMa(2)=+MaxYLIMa;
end
set(gca,'YLim',YLIMa);
hold on
plot([T_ei T_ei],YLIMa,'g--','LineWidth',1.5);
plot([T_ef T_ef],YLIMa,'g--','LineWidth',1.5);
hold off
YTICKa=get(gca,'YTick');
if( length(YTICKa)==3 )
    YTICKa=YLIMa(1):0.25*(YLIMa(2)-YLIMa(1)):YLIMa(2);
end
set(gca,'YTick',YTICKa);
grid on
leg=legend(['(' CH2 ')']);
set(leg,'FontSize',12);

%subplot(4,1,3:4)
%plot(f_zp,abs(TR),'r-','LineWidth',0.5);
%loglog(f_zp,abs(TR),'r-','LineWidth',0.5);

subplot(4,1,3:4)

% 1. 使用 'semilogy' (Y軸log，X軸linear)
semilogy(f_zp, abs(TR), 'r-', 'LineWidth', 0.5);

% 2. (請註解掉或刪除 'loglog' 和 'plot' 這一行)
% loglog(f_zp, abs(TR), 'r-', 'LineWidth', 0.5);
% plot(f_zp, abs(TR), 'r-', 'LineWidth', 0.5);

% 3. 設定 X 軸的顯示範圍
xlim([0.2 12]);  % 顯示 0.2 Hz 到 12 Hz
if (itp==2)
hold on
plot(f_zp,abs(TR_sm),'g-','LineWidth',2);
hold off
end

xlabel('Frequency (Hz)','FontSize',12);
clear YLABEL
YLABEL=['Spectral Ratio'];
HY=ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
%set(HY,'VerticalAlignment','Bottom');
% if( exist('F_tick','var') )
%     set(gca,'XLim',[F_tick(1) F_tick(end)],'XTick',F_tick);
% end
grid on

YLIM=get(gca,'YLim');
if( YLIM(1) < 0. )
    YLIM(1)=0.;
end

% YLIM(2)=10;  %% setting upper bound of TF axis

set(gca,'YLim',YLIM);

if (itp==2)
%  if (f_pre<20.)
%   text(f_pre,TR_sm_pre*1.2,['\leftarrow (' num2str(f_pre,'%5.3f') ', ' num2str(TR_sm_pre,'%5.3f') ')'],'HorizontalAlignment','left','Color','b','FontWeight','bold' ,'FontSize',12)
%  elseif (f_pre>20.)&&(f_pre<30.)
  text(f_pre,TR_sm_pre*1.05,['(' num2str(f_pre,'%5.3f') ', ' num2str(TR_sm_pre,'%5.3f') ')\rightarrow '],'HorizontalAlignment','right','Color','b','FontWeight','bold' ,'FontSize',12)
%  else
%   text(30.,TR_sm_pre,['(' num2str(f_pre,'%5.3f') ', ' num2str(TR_sm_pre,'%5.3f') ')\rightarrow '],'HorizontalAlignment','right','Color','b','FontWeight','bold' ,'FontSize',12)
%  end
end

%legend('Original','Filtered',1)
%legend('FFT',1)
if (itp ==2)
  leg=legend(['(' CH1 ') w.r.t (' CH2 ')'],'Smoothed');
else 
   leg=legend(['(' CH1 ') w.r.t (' CH2 ')']);
end
  set(leg,'FontSize',12);

% if (f_pre<=F_max*0.8)
  set(leg,'FontSize',12, 'Location', 'NorthEast');
% else
%   set(leg,'FontSize',12, 'Location', 'NorthWest');    
% end

set(gcf,'Position',[200,40,960,800]);

if (icount<10)
    CON=['00' num2str(icount)];
elseif (icount<100)
    CON=['0' num2str(icount)];
else
    CON=num2str(icount);
end


%================================================
%% output TR value
if (itp ~= 1)
 Filename=[Output_Path  'TR_(' CH1 ')wrt(' CH2 ')_' ps '.out'];
 fid1=fopen(Filename,'w');
 fprintf(fid1,'\nf       TR(f)          Phase(rad)\n'); 
 for i=1:length(f_zp)
  fprintf(fid1,'%5.3f %14.6e %10.6f\n',f_zp(i),abs(TR(i)),angle(TR(i))); 
 end
fclose(fid1);
end

if (itp == 2)
 Filename=[Output_Path 'TR[smoothed]_(' CH1 ')wrt(' CH2 ')_' ps '.out'];
 fid2=fopen(Filename,'w');
 fprintf(fid2,'\nf       TR(f)          Phase(rad)\n'); 
 for i=1:length(f_zp)
  fprintf(fid2,'%5.3f %14.6e %10.6f\n',f_zp(i),abs(TR_sm(i)),angle(TR_sm(i))); 
 end
fclose(fid2); 
end



if (itp ~= 1)
%% output TR fig
fig_out=[Output_Path 'TR_(' CH1 ')wrt(' CH2 ')_' ps '.png'];
saveas(gcf,fig_out,'png');
end

 
%% output f_pre value
% if (itp==2)
%   fprintf(fid_pre,'%9.2f %9.3f %9.3f\n',T_ei+NN*dt/2.,f_pre,TR_sm_pre); 
% end
%============================================

if (itp==1)
  icount=icount+1;
  T_ei=T_ei+NN*(1.-overlap)*dt;
else
   T_ei=T_f;
end

%%+++++++++++++++++++++
end
%%+++++++++++++++++++++

% if (itp==1)
%   fclose(fid_pre);
% end


%============================================
if (itp==1)
%============================================
TR_ave=TR_ave/(icount-1);
% TR_ave=smooth(TR_ave,1/df_zp,'moving');

[TR_ave_max,maxf_ave_ID]=max(abs(TR_ave(1:round(15/df_zp))));
f_ave_pre=f_zp(maxf_ave_ID);
TR_ave_pre=abs(TR_ave_max);

subplot(4,1,1)
plot(t,R1_fil,'b-','LineWidth',0.2);
clear TITLE
TITLE{1}=['Remark: ' ps ];
TITLE{2}=['Time Window:[ ' num2str(T_i,'%5.2f') ', ' num2str(T_f,'%5.2f') '] sec (ave) ,  T_eff=' num2str(T_effect) ' sec,  Overlap=' num2str(NN*overlap*dt) ' sec,  Filter: [' num2str(F_min) ',' num2str(F_max) '] Hz'];
title(TITLE,'FontSize',12,'Interpreter','none');
xlabel(['Time (' T_unit ')'],'FontSize',12);
clear YLABEL
YLABEL=['v(t) (' R_unit ')'];
ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
if( exist('T_tick','var') )
    set(gca,'XLim',[T_tick(1) T_tick(end)],'XTick',T_tick);
end
%axis([T_tick(1) T_tick(end) -0.005 0.005])
YLIMa=get(gca,'YLim');
if( YLIMa(1) ~= -YLIMa(2) )
  MaxYLIMa=max(abs(YLIMa));
  YLIMa(1)=-MaxYLIMa;
  YLIMa(2)=+MaxYLIMa;
end
set(gca,'YLim',YLIMa);
YTICKa=get(gca,'YTick');
if( length(YTICKa)==3 )
    YTICKa=YLIMa(1):0.25*(YLIMa(2)-YLIMa(1)):YLIMa(2);
end
set(gca,'YTick',YTICKa);
grid on
leg=legend(['(' CH1 ')']);
set(leg,'FontSize',12);


subplot(4,1,2)
plot(t,R2_fil,'b-','LineWidth',0.2);
xlabel(['Time (' T_unit ')'],'FontSize',12);
clear YLABEL
YLABEL=['v(t) (' R_unit ')'];
ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
if( exist('T_tick','var') )
    set(gca,'XLim',[T_tick(1) T_tick(end)],'XTick',T_tick);
end
%axis([T_tick(1) T_tick(end) -0.005 0.005])
YLIMa=get(gca,'YLim');
if( YLIMa(1) ~= -YLIMa(2) )
  MaxYLIMa=max(abs(YLIMa));
  YLIMa(1)=-MaxYLIMa;
  YLIMa(2)=+MaxYLIMa;
end
set(gca,'YLim',YLIMa);
YTICKa=get(gca,'YTick');
if( length(YTICKa)==3 )
    YTICKa=YLIMa(1):0.25*(YLIMa(2)-YLIMa(1)):YLIMa(2);
end
set(gca,'YTick',YTICKa);
grid on
leg=legend(['(' CH2 ')']);
set(leg,'FontSize',12);

subplot(4,1,3:4)
HL=plot(f_zp,TR_ave,'r-','LineWidth',0.5);
xlabel('Frequency (Hz)','FontSize',12);
clear YLABEL
YLABEL=['Spectral Ratio'];
HY=ylabel(YLABEL,'FontSize',12,'VerticalAlignment','Bottom');
%set(HY,'VerticalAlignment','Bottom');
if( exist('F_tick','var') )
    set(gca,'XLim',[F_tick(1) F_tick(end)],'XTick',F_tick);
end
grid on

YLIM=get(gca,'YLim');
if( YLIM(1) < 0. )
    YLIM(1)=0.;
end

set(gca,'YLim',YLIM);


% if (f_ave_pre<=15.)
%   text(f_ave_pre,TR_ave_pre,['\leftarrow (' num2str(f_ave_pre,'%5.2f') ', ' num2str(TR_ave_pre,'%5.2f') ')'],'HorizontalAlignment','left','Color','b','FontWeight','bold' ,'FontSize',12)
% elseif (f_ave_pre>15.)&&(f_ave_pre<20.)
  text(f_ave_pre,TR_ave_pre,['(' num2str(f_ave_pre,'%5.2f') ', ' num2str(TR_ave_pre,'%5.2f') ')\rightarrow '],'HorizontalAlignment','right','Color','b','FontWeight','bold' ,'FontSize',12)
% else
%   text(20.,TR_ave_pre,['(' num2str(f_ave_pre,'%5.2f') ', ' num2str(TR_ave_pre,'%5.2f') ')\rightarrow '],'HorizontalAlignment','right','Color','b','FontWeight','bold' ,'FontSize',12)
% end
    
leg=legend(['(' CH1 ') w.r.t (' CH2 ')']);
% if (f_ave_pre<=F_max*0.8)
  set(leg,'FontSize',12, 'Location', 'NorthEast');
% else
%   set(leg,'FontSize',12, 'Location', 'NorthWest');    
% end

set(gcf,'Position',[200,40,960,800]);

fig_out=[Output_Path 'TR_(' CH1 ')wrt(' CH2 ')_ave_' ps '.png'];
saveas(gcf,fig_out,'png');
% 
%% output TR_ave
Filename=[Output_Path 'TR_(' CH1 ')wrt(' CH2 ')_ave_' ps '.out'];
fid3=fopen(Filename,'w');
fprintf(fid3,'f_ave     TR_ave\n'); 
for i=1:length(f_zp)
fprintf(fid3,'%5.3f %14.6e\n',f_zp(i),TR_ave(i)); 
end

fclose(fid3);

%============================================
end
%============================================
