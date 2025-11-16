function RnC_TRF(fname1,fname2, CH1, CH2, ps,overlap,itp)
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

global F_min F_max  T_i T_f T_ei T_ef NN
global T_tick F_tick 
global Data_Path Output_Path

F_name1=[Data_Path fname1];
  [fid1,message]=fopen(F_name1,'r');
  if fid1==-1
    disp(message)
    return;
  end
  
 mult1=fscanf(fid1,'%g',1);
 dt1=fscanf(fid1,'%g',1);
 [R1,count]=fscanf(fid1,'%g',[1,inf]);
 R1=R1*mult1;
 fclose(fid1);
 nt1=count;

F_name2=[Data_Path fname2];
  [fid2,message]=fopen(F_name2,'r');
  if fid2==-1
    disp(message)
    return;
  end
  
 mult2=fscanf(fid2,'%g',1);
 dt2=fscanf(fid2,'%g',1);
 [R2,count]=fscanf(fid2,'%g',[1,inf]);
 R2=R2*mult2;
 fclose(fid2);
 nt2=count;

 if ((nt1 ~= nt2) || (dt1 ~= dt2))
    disp('R1 & R2 not compitable')
    return;
 end
 
 if (T_f > nt1*dt1)
      T_f= nt1*dt1;
 end  

XQ_TRF(R1,R2,nt1,dt1,CH1,CH2,ps,overlap,itp);
  
end