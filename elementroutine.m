function[k,F,stress,overstress,strain]=elementroutine(r1,r2,c,N,zi,u,T,Q,dt,delta_u,prv_overstress)
%%calculating strain displacemnet matrix matrix
B=[-1/(r2-r1),1/(r2-r1);(1-zi)/(r1*(1-zi)+r2*(1+zi)),(1+zi)/(r1*(1-zi)+r2*(1+zi))];
%% Jacobian
J=(r2-r1)/2;
r=[r1;r2]; 
strain=B*u;
stress_linear=c*strain;                %% stress for the linear case
[c,overstress]=materialroutine(T,c,Q,B,dt,delta_u,prv_overstress);
stress=stress_linear+overstress;       %% stress for the non linear case
k=2*B'*c*B*(N'*r)*(J);                 %% element stiffness matrix
F=2*B'*stress*(N'*r)*(J);              %% element internal forces
end 
