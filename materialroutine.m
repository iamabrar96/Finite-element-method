 function [ct,overstress]=materialroutine(T,c,Q,B,dt,delta_u,prv_overstress)
delstrain=B*delta_u;           
vol_strain=(sum(delstrain)/3)*[1;1];
dev_strain=delstrain-vol_strain;          %% calculation of deviotoric strain 
%% internal stack variable(overstress)
overstress= (1/(1+(dt/(2*T)))) * ((prv_overstress*(1-(dt/(2*T)))) + Q*dev_strain);
%% material tangent stiffness matrix
ct=c+(1/(1+(dt/(2*T))))*Q*[2/3,1/3;-1/3,2/3]; 
 end  