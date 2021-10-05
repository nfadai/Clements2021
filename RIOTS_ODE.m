function RIOTS_ODE(Rec,Def,r0,b0,tend)
%% Spatially-uniform continuum limit considered in Section 4.2.1
% Use this code and RIOTS_ABM_spatiallyuniform.m to reproduce Figure 11

% kr = recruitment rate
% kd = defection rate
% r0 = initial rioter density
% b0 = initial bystander density
% tend = total run time

%% Growth rate equations

function dudt=f(t,u)
    dr=u(2)*(Rec(1)*(1-u(1)).^4 + Rec(2)*4*u(1)*(1-u(1)).^3 + Rec(3)*6*u(1).^2*(1-u(1)).^2+...
         Rec(4)*4*u(1).^3*(1-u(1))+Rec(5)*u(1).^4) - ...
         u(1)*(Def(1)*(1-u(2)).^4 + Def(2)*4*u(2)*(1-u(2)).^3 + Def(3)*6*u(2).^2*(1-u(2)).^2+...
         Def(4)*4*u(2).^3*(1-u(2))+Def(5)*u(2).^4);
   
    db=-dr;
dudt=[dr;db];
end

%% Initial Conditions

u0=[r0;b0];

%% Simulation time run

trun=[0 tend];

%% Solving the ODE system

[t,u]=ode45(@f,trun,u0);

%% Plotting trajectories

figure(501)
plot(t,u(:,1),'--','linewidth',2,'color',[0.64,0.08,0.18])
hold on
plot(t,u(:,2),'c--','linewidth',2)
%hold off
xlabel('t')
ylabel('r(t), b(t)')
legend({'r(t)' 'b(t)'})
set(gca,'FontSize',20)
end