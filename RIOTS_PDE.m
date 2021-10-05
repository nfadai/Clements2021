function RIOTS_PDE(D,Rec,Def,X,tend,m0r,m0b)
%% Growth rate Equations
% Writes the growth rates for the continuum limit system
% u(1) is (mean) rioter occupancy
% u(2) is (mean) bystander occupancy

m=0;

    function [c,f,s]=pdefun(x,t,u,dudx)
        c=[1;1];
        f=D*[dudx(1).*(1-u(2))+u(1).*dudx(2);dudx(2).*(1-u(1))+u(2).*dudx(1)];
            dr=u(2)*(Rec(1)*(1-u(1)).^4 + Rec(2)*4*u(1)*(1-u(1)).^3 + Rec(3)*6*u(1).^2*(1-u(1)).^2+...
         Rec(4)*4*u(1).^3*(1-u(1))+Rec(5)*u(1).^4) - ...
         u(1)*(Def(1)*(1-u(2)).^4 + Def(2)*4*u(2)*(1-u(2)).^3 + Def(3)*6*u(2).^2*(1-u(2)).^2+...
         Def(4)*4*u(2).^3*(1-u(2))+Def(5)*u(2).^4);
   
    db=-dr;
s=[dr;db];
    end

%% Initial Conditions
% Equivalent to initial conditions of ABM with rioters and bystanders at
% opposite ends of the street (street) and rioters in a central strip (central)

% Street initial conditions

%    function u0=icfun(x)
%        if x>=80 && x<=90
%            u0=[1;0];
%        elseif x>=110 && x<=120
%            u0=[0;1];
%        else
%            u0=[0;0];
%        end
%    end

% Central initial conditions

    function u0=icfun(x)
        if x>=60 && x<=80||x>=120 && x<=140
            u0=[0;m0b];
        elseif x>=90 && x<=110
            u0=[m0r;0];
        else
            u0=[0;0];
        end
    end

%% Boundary Conditions
% Neumann (zero flux) boundary conditions

    function [pl,ql,pr,qr]=bcfun(xl,ul,xr,ur,t)
        pl=[0;0];
        ql=[1;1];
        pr=[0;0];
        qr=[1;1];
    end

%% Spatial and time domains
% Defines 1d x-space 0<=x<=200 and time to run until 'tend'

xval=linspace(1/2,X-1/2,1000);
tval=linspace(0,tend,1000);

%% Solve the PDE
% Solves the continuum limit equations

sol=pdepe(m,@pdefun,@icfun,@bcfun,xval,tval);

r=sol(:,:,1); % Captures rioter distribution solution
b=sol(:,:,2); % Captures bystander distribution solution

%% Plot solutions in x-space

figure(86)
plot(xval,r(end,:),'--','linewidth',2,'color',[0.64,0.08,0.18])
hold on
plot(xval,b(end,:),'c--','linewidth',2)
xlim([0,X])
ylim([0,1])
xlabel('x')
ylabel('Agent Densities')
legend({'r(x,t)' 'b(x,t)'})
hold off
set(gca,'Fontsize',20)
end
