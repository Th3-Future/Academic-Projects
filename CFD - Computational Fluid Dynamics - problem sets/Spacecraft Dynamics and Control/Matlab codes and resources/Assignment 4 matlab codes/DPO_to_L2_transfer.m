clear
clc
load A4orbits.mat
mu=0.01213;         %mass ratio for moon

coeff=[-1 mu+2 -2*mu-1 mu+1 2*mu-2 1-mu];       %coefficients for polynomial to solve for L2 location
C1=roots(coeff);
D=imag(C1);
D=abs(D);
[~,I]=min(D);

rs1_c=C1(I);
xL2_c=rs1_c-mu;
rs2_c=rs1_c-1;
yL2_c=0;

%giving initial conditions
xL2=X0L2(1,1);
yL2=X0L2(2,1);
dxh=X0L2(3,1);
dyh=X0L2(4,1);
epsilon=2.6e-4;

%Jacobi Constant
rs1h=sqrt((xL2+mu)^2+yL2^2);
rs2h=sqrt((xL2-1+mu)^2+yL2^2);
C=2*((xL2^2+yL2^2)/2+(1-mu)/rs1h+mu/rs2h)-dyh^2-dxh^2

%% Full state vector at L2, equal to the given initial condition
XL2=X0L2;
options=odeset('RelTol',1e-12,'AbsTol',1e-15);

%%compute eigenvalues and eigenvectors of the monodromy matrix 
%%for planar dynamics of L2
X0=[XL2;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
[Tsp,Xsp] = ode45(@pcrtbp_dyn_stm,[0 TpL2],X0,options);
%extract state transition matrix at end point
[ros,~]=size(Tsp);
%monodromy matrix
STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).']; 
%eigen decomposition of monodromy matrix
[V,D]=eig(STM)

%initialize array for the poincare map
y1=zeros(ros,1); 
y1dot=zeros(ros,1);

%figure definition
fig1 = figure(1) ;ax1 = axes ;hold(ax1,'on')
fig2 = figure(2) ;ax2 = axes ; hold(ax2,'on')
fig3 = figure(3) ;ax3 = axes ; hold(ax3,'on')

plot(ax1,1-mu,0,'ko')
plot(ax1,xL2_c,0,'r+')

%For each state on the L2 Orbit, plot the stable manifold
for test=1:1:ros
    stm=[Xsp(test,5:8).' Xsp(test,9:12).' Xsp(test,13:16).' Xsp(test,17:20).'];
    u_s=real(V(:,2));   %eigenvector for stable manifold
    ratio=(stm*u_s)/norm(stm*u_s);
    Xup0=Xsp(test,1:4).'+epsilon*(ratio);
    options1=odeset('RelTol',1e-12,'AbsTol',1e-15,'Events',@myEvent);
    [Tup,Xup,teL,xe] = ode45(@pcrtbp_dyn2,[0 TpL2],Xup0,options1);
    y1(test)=xe(1,2);
    y1dot(test)=xe(1,4);
    plot(ax1,Xup(:,1),Xup(:,2),'--b','LineWidth',0.001)
end

%% Full state vector at DPO, equal to the given initial condition
Xm=X0m;

%%compute eigenvalues and eigenvectors of the monodromy matrix 
%%for planar dynamics of moon DPO
X0=[Xm;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
[Tsp,Xspm] = ode45(@pcrtbp_dyn_stm,[0 Tpm],X0,options);
%extract state transition matrix at end point
[ros,~]=size(Tsp);
Xf=Xspm(ros,1:4).';
%monodromy matrix
STM=[Xspm(ros,5:8).' Xspm(ros,9:12).' Xspm(ros,13:16).' Xspm(ros,17:20).']; 
%eigen decomposition of monodromy matrix
[V,D]=eig(STM)

%initialize array for the poincare map
y2=zeros(ros,1); 
y2dot=zeros(ros,1);

% For each state on the DPO Orbit, plot the unstable manifold
for test=1:1:ros
    stm=[Xspm(test,5:8).' Xspm(test,9:12).' Xspm(test,13:16).' Xspm(test,17:20).'];
    u_u=real(V(:,1));   %eigenvector for unstable manifold
    ratio=(stm*u_u)/norm(stm*u_u);
    Xup0=Xspm(test,1:4).'+epsilon*(ratio);
    options1=odeset('RelTol',1e-12,'AbsTol',1e-15,'Events',@myEvent);
    [Tup,Xup,teM,xe] = ode45(@pcrtbp_dyn,[0 10],Xup0,options1);
    plot(ax1,Xup(:,1),Xup(:,2),'--r','LineWidth',0.001)
    y2(test)=xe(1,2);
    y2dot(test)=xe(1,4);
end

%%
%Poincare Map plot
plot(ax2,y1dot(:),y1(:),'-b')
plot(ax2,y2dot(:),y2(:),'-r')

%% Free Transfer from DPO to L2
%Finds all intersection, but we have only one true intersection point
[dyi,yi] = polyxpoly(y1dot(:),y1(:),y2dot(:),y2(:));
scatter(ax2,dyi(5),yi(5),'*k') %plots the point

%calculate x'
x=xe(1);
y=yi(5);
dy=dyi(5);
rs1=sqrt((x+mu)^2+y^2);
rs2=sqrt((x-1+mu)^2+y^2);
dx=sqrt(2*(((x^2+y^2)/2)+((1-mu)/rs1)+mu/rs2)-dy^2-C); %from note, eqn(78)
%Found initial condition for transfer
x0=[x;y;dx;dy];
%plot L2 and DPO orbits
plot(ax3,Xsp(:,1),Xsp(:,2),'-b')
plot(ax3,Xspm(:,1),Xspm(:,2),'-r')

%%compute transfer
[TspL,XspL]=ode45(@pcrtbp_dyn,[0 teL],x0,options);%forward in time
[TspM,XspM]=ode45(@pcrtbp_dyn2,[0 teM],x0,options);%backward in time
%plot transfers
plot(ax3,XspL(:,1),XspL(:,2),'--b')
plot(ax3,XspM(:,1),XspM(:,2),'--r')
save DPO_to_L2.mat TspM TspL XspM XspL;

%% Show all Plots
%%(N/B:holding all plots till here, reduces simulation time)

figure(1)
leg1=legend('$Moon$','$L_2$','L2 stable Manifold','DPO unstable Manifold');
set(leg1,'Interpreter','latex');
title('Stable & Unstable Manifolds for Periodic L2 and LunarDPO Orbits');
xlabel('x');
ylabel('y');
hold off

figure(2)
leg1=legend('$L_2$ stable','Lunar DPO unstable');
set(leg1,'Interpreter','latex');
title('Poincare Maps');
xlabel('ydot');
ylabel('y');
hold off

figure(3)
leg1=legend('$L_2$ Orbit','Lunar DPO Orbit');
set(leg1,'Interpreter','latex');
title('Transfer');
xlabel('x');
ylabel('y');
hold off

%N/B for plot colour, we are performing tranfer from the red to the blue one.