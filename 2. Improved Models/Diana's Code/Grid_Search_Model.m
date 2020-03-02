function [SI,d2]=Grid_Search_Model(d2,val,Treal,Greal,Ireal,timestep,Dinput,Iinput)

format long

val_temp=val;

jmax=length(d2);

SI=zeros(jmax,1);
errg=zeros(jmax,1);

for j=1:jmax
    val_temp.d2=d2(j);
        
    [SI_calc]=SI_Model(val_temp,Treal,Greal,Ireal,timestep,Dinput,Iinput);
    SI(j,1)=SI_calc;
    val_temp.SI=SI_calc;
        
    [output,T_GIQ]=GIQ_Model(val_temp,Greal,Ireal,timestep,Dinput,Iinput);

    G=interp1(T_GIQ,output(:,5),Treal);
    errg(j,1)=mean(((abs(Greal-G))./Greal)*100);
end
minimum=min(min(errg));
[j]=find(errg==minimum);
SI=SI(j,1);
d2=d2(j);

end


% function [SI,da,d2]=Grid_Search_Model(da,d2,val,Treal,Greal,Ireal,timestep,Dinput,Iinput)
% 
% format long
% 
% val_temp=val;
% 
% imax=length(da);
% jmax=length(d2);
% 
% SI=zeros(imax,jmax);
% errg=zeros(imax,jmax);
% 
% for i=1:imax
%     for j=1:jmax
%         val_temp.da=da(i);
%         val_temp.d2=d2(j);
%         
%         [SI_calc]=SI_Model_delay(val_temp,Treal,Greal,Ireal,timestep,Dinput,Iinput);
%         SI(i,j)=SI_calc;
%         val_temp.SI=SI_calc;
%         
%         [output,T_GIQ]=GIQ_Model_delay(val_temp,Greal,Ireal,timestep,Dinput,Iinput);
% 
%         G=interp1(T_GIQ,output(:,6),Treal);
%         errg(i,j)=mean(((abs(Greal-G))./Greal)*100);
%     end
% end
% minimum=min(min(errg));
% [i,j]=find(errg==minimum);
% SI=SI(i,j);
% da=da(i);
% d2=d2(j);
% 
% end