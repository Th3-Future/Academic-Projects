clear
clc

mu=0.01213;         %mass ratio for moon
coeff=[-1 mu+2 -2*mu-1 mu+1 2*mu-2 1-mu];       %coefficients for polynomial to solve for L2 location
C1=roots(coeff);
D=imag(C1);
D=abs(D);
[~,I]=min(D);
rs1_c=C1(I);
xL2_c=rs1_c-mu;

%% Load transfer result variables
load DPO_to_L2.mat;
load L2_to_DPO.mat;
%% creating n-number of segment from the transfer states
rosL=1746;   % index close to surface of section crossing (x=1.075)
rosM=1787;
%four segments
s1=flip(XspM(1:rosM,1:2));
s2=XspL(1:rosL,1:2);
s3=flip(XspL2(1:rosL,1:2));
s4=XspM2(1:rosM,1:2);
%and corresponding timespan
t1=TspM(1:rosM);
t2=TspL(1:rosL)+t1(end);
t3=TspL2(1:rosL)+t2(end);
t4=TspM2(1:rosM)+t3(end);
%and corresponding velocity
v1=flip(XspM(1:rosM,3:4));
v2=XspL(1:rosL,3:4);
v3=flip(XspL2(1:rosL,3:4));
v4=XspM2(1:rosM,3:4);

%selecting patch point for each segments
s1(1,2)=0; s4(end,2)=0;%setting point 1 and 5 on x-axis crossing
pp12=s1(1,:);
pp23=s2(1,:);
pp34=s3(1,:);
pp45=[s4(1,:);s4(end,:)];
%timespane for the each segment
t12=t1(1,:);
t23=t2(1,:);
t34=t3(1,:);
t45=[t4(1,:);t4(end,:)];
%velocity for the each segment
v12=v1(1,:);
v23=v2(1,:);
v34=v3(1,:);
v45=[v4(1,:);v4(end,:)];

PP=[pp12,t12;pp23,t23;pp34,t34;pp45,t45].';
Vg=[v12;v23;v34;v45].';
%% Multishooting method
[X0,t0,tf,PPf,VV,YY]=pcrtbpmultshoot(PP,Vg);

%%compute the dynamics of obtained initial state
options=odeset('RelTol',1e-12,'AbsTol',1e-15);
[Tf,Xf]=ode45(@pcrtbp_dyn,[t0 tf],X0,options);

%plot new orbit about DPO and L2
fig1 = figure(1) ;ax1 = axes ;hold(ax1,'on')
plot(ax1,1-mu,0,'mo')
plot(ax1,xL2_c,0,'g+')
plot(ax1,Xf(:,1),Xf(:,2),'-r','LineWidth',1)
leg1=legend('$Moon$','$L_2$');
set(leg1,'Interpreter','latex');
title('New Orbit about Lunar DPO and L2');
xlabel('x');
ylabel('y');
hold off
