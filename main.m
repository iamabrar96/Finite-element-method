%% Abrar Hyder Mohammed
%% 65092
%% NLFEM Assignemnet 2020
clc;
clear all;
prompt = '[1]:Non-linear [0]:linear?: ';
x = input(prompt);
if x==0
    Q=0;
else
    Q=35000;
end
%% Input parameters:

r1=40;                                   %% inner radius
r2=80;                                   %% outer raius
E=70000;                                %% youngs modulus              
T=2;                                     %% Time scale  in sec
p=50;                                    %% maximum pressure
tl=4;                                    %% initial load time
tf=20;                                   %% final load time
v=0.25;                                  %% poison's ratio
nelem=10;                                %% number of elements
dt=(1/nelem);
alpha=((E)/((1+v)*(1-2*v)));
beta=[1-v,v;v,1-v];
c=alpha*beta;                            %% linear elasticity
zi=0;                                    %% quadrature with one gauss point
n1=(1-zi)/2;                             %% shape functions 
n2=(1+zi)/2;
N=[n1;n2];
kg=zeros(nelem+1,nelem+1);               %% initialising global stiffness matrix
rnodes=meshGenerator(r1,r2,nelem);       %% all the required number of nodes are generated 
u=zeros(nelem+1,1);                      %% initialising global nodal dispplacement vector
delta_u = zeros(nelem+1,1);              %% initialising updated displacemnet 
overstress = zeros(nelem,1);                       
strain_g = zeros(nelem,2);
stress_g= zeros(nelem,2);
prv_overstress=zeros(nelem,1);
time_history=[];                      %% storing the time_history values in the array for plotting
displace_history=[];
for t=dt:dt:tf
    if t<=tl
        fext=zeros(nelem+1,1);
        ld_scaling = (1/tl)*t;
        fext(1) = p*r1*ld_scaling;       %% load scaling operation
        fext(1)=fext(1);
    else
        fext=zeros(nelem+1,1);
        fext(1)=p*r1;
        fext(1)=fext(1);
    end
       
    iter=0;
    while true
        fint=zeros(nelem+1,1);
        fext=zeros(nelem+1,1);
        fext(1) = p*r1*ld_scaling;    
        fext(1)=fext(1);
        kg=zeros(nelem+1,nelem+1);
       
        for i=1:length(rnodes)-1
            %% calculating element stiffness and internal forces from element routine
            [k,F,stress,overstress,strain]=elementroutine(rnodes(i),rnodes(i+1),c,N,zi,u(i:i+1),T,Q,dt,delta_u(i:i+1),prv_overstress(i));
            prv_overstress(i)=overstress(1);
            stress_g(i,1)= stress(1);
            stress_g(i,2) = stress(2);
            %% Assembly of global stiffness matrix and internal forces 
            kg(i:i+1,i:i+1)=kg(i:i+1,i:i+1)+k;    
            fint(i:i+1)=fint(i:i+1)+F;             
            
        end
        %% solving for system of linear equations
        delta_u = linsolve(kg,fext-fint);   
        u = u+delta_u;
        if (max(abs(fint-fext)))<0.005*(max(abs(fint)))|| (max(abs(delta_u)))<0.005*(max(abs(u)))
            disp(iter)
            break;         
        else  
            iter=iter+1;
        end   
    end
   time_history=[time_history,t];               %% appending time_history with time 
   displace_history=[displace_history,u];
end
                      %% Plotting
                      
if x==0
    figure(1);
    Analyticalsol=Analytical(v,p,E,r1,r2,rnodes);
    plot(rnodes,Analyticalsol,'r--o',rnodes,u,'b-*');
    title('Analytical solution vs Numerical solution')
    xlabel(' r (mm)');
    ylabel('u_{r} (mm)');
    legend('Analytical','Numerical');
    legend('boxoff');
    savefig('linear.fig');
    print('Linear','-dpng')

    figure(2);
    plot(time_history,displace_history(end,:),'b-','linewidth',2);
    title('Time history of widening of the pipe','Interpreter', 'tex')
    xlabel('Time t (s)');
    ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
    legend('last node @ elastic case','location','southeast');
    savefig('time_history1.fig');
    print('time_history1','-dpng')

else
    figure(3);
    Analyticalsol=Analytical(v,p,E,r1,r2,rnodes);
    plot(rnodes,Analyticalsol,'r--o',rnodes,u,'b-*');
    title('Analytical solution vs Numerical solution')
    xlabel(' r (mm)');
    ylabel('u_{r} (mm)');
    legend('Analytical','Numerical');
    legend('boxoff');
    savefig('nonlinear.fig');
    print('nonlinear','-dpng')

    figure(4);
    hold on;
    plot(time_history,displace_history(end,:),'r-','linewidth',2);
    hold off;
    title('Time history of widening of the pipe','Interpreter', 'tex')
    xlabel('Time t (s)');
    ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
    legend('last node @ viscoelastic case','location','SouthEast');
    legend('boxoff'); 
    savefig('time_history2.fig');
    print('time_history2','-dpng')

    figure(5);
    subplot(2,1,1);
            plot(rnodes(1:end-1),stress_g(:,1),'b-*','linewidth',1);
            title('{\sigma}_{r r} vs  r','Interpreter', 'tex');
            xlabel('radius r (mm)');
            ylabel('{\sigma}_{r r}(MPa)','Interpreter', 'tex');
            legend('stress rr','location','southeast');
        subplot(2,1,2);
            plot(rnodes(1:end-1),stress_g(:,2),'g-o','linewidth',1);
            title('{\sigma}_{\phi \phi} vs  r','Interpreter', 'tex');
            xlabel('radius r (mm)');
            ylabel('{\sigma}_{\phi \phi} (MPa)','Interpreter', 'tex');
            legend('stress ff');
   savefig('stress rr&ff.fig');
   print('stress rr&ff','-dpng')

    
end




