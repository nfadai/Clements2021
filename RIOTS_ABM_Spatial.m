function RIOTS_ABM_Spatial(Ttot,P)
%% Runs P simulations of ABM for riots model with recruitment and defection
% kmr = movement rate for rioters (unbiased in direction)
% kmb = movement rate for bystanders (unbiased in direction)
% kr = recruitment rate
% kd = defection rate
% Ttot = Overall simulation time
% P = Number of identically-prepared simulations to run
%% Setup and Initial Conditions

%rng(5)
X=200; % Number of x-ordinates
Y=20; % Number of y-ordinates
set(0, 'DefaultLineLineWidth', 2);
set(gca,'FontSize',20)
if nargin==0
    Ttot=8;
    P=20;
end
m0r=0.5; % Initial rioter density (for spatial-uniformity)
m0b=0.5; % Initial bystander density (for spatial-uniformity)

K=0;

AA=0;BB=0;CC=0;DD=0;
Rec=[0 1 1 1 1];
% DEFINE INDIVIDUAL RECRUITMENT RATES HERE
Def=[0 1 1 1 1];
    (CC+3*K*DD)*[0 0 0 1/4 1] - DD*[0 0 0 0 1]; % DEFINE INDIVIDUAL DEFECTION RATES HERE



%% Nearest neighbour index structure

XX=zeros(Y,X); % To contain row index i for each (i,j)
YY=zeros(Y,X); % To contain column index j for each (i,j)


if min(min(Rec),min(Def))<0
    error('ERROR: negative individual-level rates chosen')
    
end
Pmr=100*max(max(Rec),max(Def)); % Moving prob/rate for rioters
Pmb=100*max(max(Rec),max(Def)); % Moving prob/rate for bystanders

for i=1:Y
    for j=1:X
        XX(i,j)=i;
        YY(i,j)=j;
    end
end
S=(X*Y+1)*ones(X*Y,4); % Set up matrix for neighbours, (XY+1) where less neighbours

for i=1:Y
    for j=1:X
        Z=(XX-XX(i,j)).^2+(YY-YY(i,j)).^2-1; % Metric for distance from (i,j)
        V=find(Z<=1e-1 & Z>-1); % Chooses sites neighbouring (i,j)
        S(i+(j-1)*Y,1:length(V))=V; % Lists neighbour indexes for (i,j) indexed by J, along row J
    end
end

S(S==0)=X*Y+1; % Buffer

parfor p=1:P
    CCR=zeros(Y,X); % Empty lattice for tracking rioters (0s will denote empty sites)
    CCB=zeros(Y,X); % Empty lattice for tracking bystanders
    
    % Spatially-uniform initial conditions
    
    %CCR = double(rand(Y,X)<m0r); % Places the rioters according to density
    %II=find(CCR==0); % Finds indices of where bystanders are allowed to be potentially placed
    %CCB(II(rand(1,length(II))<m0b/(1-m0r)))=1; % Places the bystanders acc. to density (and where no rioters)
    
    
    % Street initial conditions
    
    % CCR(:,81:90)=1; % Initial conditions (no spatial-uniformity) Left end of avenue full of rioters (denoted by 1s)
    % CCB(:,111:120)=1; % Initial conditions (no spatial-uniformity) Right side of avenue full of bystanders (denoted by 1s)
    
    
    % Central initial conditions
    
    CCR(:,91:110)=double(rand(size(CCR(:,91:110)))<m0r); % Initial conditions 2 (no spatial-uniformity) Central stripe of rioters
    CCB(:,61:80)=double(rand(size(CCB(:,61:80)))<m0b); % Initial conditions 2 (no spatial-uniformity) Bystander crowd to left
    CCB(:,121:140)=double(rand(size(CCB(:,121:140)))<m0b); % Initial conditions 2 (no spatial-uniformity) Bystander crowd to right
    
    
    QRstart=0; % Counts the total number of rioters initially
    QBstart=0; % Counts the total number of bystanders initially
    for y=1:Y
        for x=1:X % Counts over every site on the lattice
            if CCR(y,x)==1 % There's a rioter in that site
                QRstart=QRstart+1; % Counts a rioter
            elseif CCB(y,x)==1 % There's a bystander in that site
                QBstart=QBstart+1; % Counts a bystander
            else
            end
        end
    end % QRstart and QBstart now have the total numbers of rioters and bystanders
    
    
    %     CRstart=CCR; % Initial profile of rioters
    %     CBstart=CCB; % Initial profile of bystanders
    %
    
    
    %% Running Simulations
    
    % Runs a simulation for each p from same initial conditions
    
    T=0; % Time starts at zero
    j=1; % Index to be used for steps
    
    tau=0; % Actually, to be a vector of cumulative timesteps
    QR=QRstart; % Actually, to be a vector of total numbers of rioters
    QB=QBstart; % Actually, to be a vector of total numbers of bystanders
    
    CR=CCR; % To manipulate the lattice
    CB=CCB;
    
    D=randi([1,4],1,10000);
    JJ=randi([1,X*Y],1,10000); % Vector (length 10000) of [1,XY] randm integers
    U1=rand(1,10000); % First vector (length 10000) of Unif(0,1)s
    U2=rand(1,10000); % Second vector (length 10000) of Unif(0,1)s
    %U3=rand(1,10000); % Third vector (length 10000) of Unif(0,1)s
    %U4=rand(1,10000); % Fourth vector (length 10000) of Unif(0,1)s
    y=1; % Arbitrary index
    %   QRend=QRstart; % To be updated after step
    %  QBend=QBstart; % To be updated after step
    
    CR(Y,X+1)=0;
    CB(Y,X+1)=0;
    %% Giant Gillespie Algorithm
    
    while T<Ttot % While simulation has time left to run
        J=JJ(y);
        while CB(J)==0 && CR(J)==0
            y=y+1;
            if y==10001
                D=randi([1,4],1,10000);
                JJ=randi([1,X*Y],1,10000);
                U1=rand(1,10000);
                U2=rand(1,10000);
                %U3=rand(1,10000);
                %U4=rand(1,10000);
                y=1;
            end
            J=JJ(y);
        end
        % First pick agent, then pick agent event later
        if CR(J)==1
            % Agent chosen will be a rioter
            Agent=1;
            %             Ind = find(CR); % Finds lattice sites occupied by rioters
            %             J=Ind(randi(length(Ind))); % Pick a rioter
            v=S(J,:); % Contains indexes of neighbours of J
            
            Nr=sum(CR(v)); % Number of rioter neighbours
            Nb=sum(CB(v)); % Number of bystander neighbours
            Pr=Rec(Nr+1); % Picks the relevent individual-level rec. rate
            Pd=Def(Nb+1); % Picks the relevant individual-level def. rate
        else
            % Agent chosen will be a bystander
            Agent=2;
            %             Ind = find(CB); % Finds lattice sites occupied by bystanders
            %             J=Ind(randi(length(Ind))); % Pick a bystander
            v=S(J,:); % Contains indexes of neighbours of J
            Nr=sum(CR(v)); % Number of rioter neighbours
            Nb=sum(CB(v)); % Number of bystander neighbours
            Pr=Rec(Nr+1); % Picks the relevent individual-level rec. rate
            Pd=Def(Nb+1); % Picks the relevant individual-level def. rate
        end
        
        prop=(Pmr+Pd)*QR(j)+(Pmb+Pr)*QB(j); % Computes (approximate) propensity of the whole system
        tau(j+1)=tau(j)-log(U1(y))/prop; % Next entry of tau is time after timestep
        
        QR(j+1)=QR(j); % Creates a new element for defining later
        QB(j+1)=QB(j); % Creates a new element for defining later
        
        if Agent==1 && U2(y)<Pmr/(Pmr+Pd)  % If rioter chosen, in the event a rioter moves
            Jnext=v(D(y)); % Chooses the direction to move with no bias
            if  Jnext<(X*Y+1) && CR(Jnext)==0 && CB(Jnext)==0 % Checks site is in matrix and not occupied
                CR(Jnext)=1; % Neighbouring site now contains rioter
                CR(J)=0; % Site originally selected now empties
            else % Movement aborted
            end
            
        elseif Agent==2 && U2(y)<Pmb/(Pmb+Pr) % If bystander chosen, in the event a bystander moves
            Jnext=v(D(y)); % Chooses the direction to move with no bias
            if  Jnext<(X*Y+1) && CR(Jnext)==0 && CB(Jnext)==0 % Checks site is in matrix and not occupied
                CB(Jnext)=1; % Neighbouring site now contains bystander
                CB(J)=0; % Site originally selected now empties
            else % Movement aborted
            end
            
        elseif Agent==2 % If bystander chosen, in the event a bystander recruits
            CB(J)=0; % Site originally selected loses bystander...
            CR(J)=1; % ... who turns into a rioter
            QR(j+1)=QR(j+1)+1; % Gain rioter
            QB(j+1)=QB(j+1)-1; % Lose bystander
            
            
        elseif Agent==1 % If rioter chosen, in the event a rioter defect
            
            CR(J)=0; % Site originally selected loses rioter...
            CB(J)=1; % ... who turns into a bystander
            QR(j+1)=QR(j+1)-1; % Lose rioter
            QB(j+1)=QB(j+1)+1; % Gain bystander
        end
        
        T=tau(j+1) ; % Updating time for the next iteration
        %        QRend=QR(j+1); % Updating number of rioters after step
        %        QBend=QB(j+1); % Updating number of bystanders after step
        
        y=y+1;
        if y==10001
            U1=rand(1,10000);
            U2=rand(1,10000);
            %U3=rand(1,10000);
            D=randi([1,4],1,10000);
            %U4=rand(1,10000);
            
            JJ=randi([1,X*Y],1,10000);
            y=1;
        end
        
        % CRM(j+1,:)=mean(CR,1);
        % CBM(j+1,:)=mean(CB,1);
        j=j+1; % Increase step counter to proceed to next step
        
        
        %         figure(40)
        %         plot3(YY(CR==1),XX(CR==1),2*CR(CR==1),'ro','markersize',5,'markerfacecolor','r')
        %         view(0,90)
        % pause(0.001)
        
        
        if CR(Y,X+1)==1 || CB(Y,X+1)==1
            error('NULL POINTER ERROR')
        end
    end
    % Agent snapshots
    figure(15)
    plot3(YY(CR(1:Y,1:X)==1),XX(CR(1:Y,1:X)==1),2*CR(CR(1:Y,1:X)==1),'ro','markersize',5,'markerfacecolor','r')
    
    hold on
    plot3(YY(CB(1:Y,1:X)==1),XX(CB(1:Y,1:X)==1),CB(CB(1:Y,1:X)==1),'bo','markersize',5,'markerfacecolor','b')
    delete(findall(findall(gcf,'Type','axe'),'Type','text'))
    hold off
    set(gca,'FontSize',20)
    box on
    view(0,90)
    xlim([0,X])
ylim([0,Y])
    pause(0.01)
    
    
    %% Interpolating and avergaing profiles
    
    
    QBi(p,:)=mean(CB(1:Y,1:X),1); % Store interpolation vector for p-th simulation
    
    QRi(p,:)=mean(CR(1:Y,1:X),1); % Store interpolation vector for p-th simulation
    
end

t=linspace(0,Ttot,1000);
QQb=mean(QBi,1); % Average over simulations
QQr=mean(QRi,1); % Average over simulations
QQbs=std(QBi,1); % Average over simulations
QQrs=std(QRi,1); % Average over simulations
figure(50)
hold off
plot(.5:X-.5,QQr,'r',.5:X-.5,QQb,'b','LineWidth',2) % Plot averaged total densities
xlim([0,X])
ylim([0,1])
xlabel('x')
ylabel('Agent Densities')
legend({'<R(x,1)>' '<B(x,1)>'})
hold off
set(gca,'Fontsize',20)

RIOTS_PDE(Pmr/4,Rec,Def,X,Ttot,m0r,m0b)

