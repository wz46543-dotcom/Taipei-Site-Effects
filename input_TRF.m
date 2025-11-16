clear global
global F_min F_max  T_i T_f T_ei T_ef NN
global T_tick F_tick 
global Data_Path Output_Path

if( ~exist('DataPath','var') )
     DataPath=input('DataPath=','s');
end
Data_Path=DataPath;
if( ~exist('OutputPath','var') )
     OutputPath=input('OutputPath=','s');
end
Output_Path=OutputPath;

if( ~exist('NP','var') )
     NP=1024;
end
NN=NP;

if( ~exist('Fmin','var') )
   F_min=0.;
  else
   F_min=Fmin;
end
%disp(['F_min=' num2str(F_min)])
if( ~exist('Fmax','var') )
%   F_max=octF(end);
    F_max=0.;
  else
   F_max=Fmax;
end

if( ~exist('Ti','var') )
   T_i=0;
else
   T_i=Ti;   
end
if( ~exist('Tf','var') )
   T_f=T_i+60;
else
   T_f=Tf;  
end
if( ~exist('Tei','var') )
  T_ei=0.;
  T_ef=0.;
else
  T_ei=Tei;
  T_ef=Tei+5.12;
end
%
%if( ~exist('Tef','var') )  
%  T_ef=0.;
%else
%  T_ef=Tef;
%end
%T_I=T_i;
%T_F=T_f;
%T_EFFECT=T_F-T_I;
clear t R Reff Ti Tf 
%=========================================================================================
if((T_f-T_i)>=1000)
     DT_tick=100;
elseif(((T_f-T_i)>=500)&(T_f-T_i)<1000)
     DT_tick=50;
elseif(((T_f-T_i)>=200)&(T_f-T_i)<500)
     DT_tick=20;
elseif(((T_f-T_i)>=100)&(T_f-T_i)<200)
     DT_tick=10;
elseif (((T_f-T_i)>=40)&((T_f-T_i)<100))
     DT_tick=5;
elseif (((T_f-T_i)>=20)&((T_f-T_i)<40))
     DT_tick=2;
elseif (((T_f-T_i)>=10)&((T_f-T_i)<20))
     DT_tick=1;
elseif (((T_f-T_i)>=5)&((T_f-T_i)<10))
     DT_tick=0.5;
elseif (((T_f-T_i)>=2)&((T_f-T_i)<5))
     DT_tick=0.2;
else
     DT_tick=0.1;  
end

if( ~exist('Ttick','var') )
   T_tick=fix(T_i):DT_tick:fix(T_f+0.5);
  else
   T_tick=Ttick;
end

if( ~exist('Ftick','var') )
  F_tick=0:2:20;
else
   F_tick=Ftick;
end

