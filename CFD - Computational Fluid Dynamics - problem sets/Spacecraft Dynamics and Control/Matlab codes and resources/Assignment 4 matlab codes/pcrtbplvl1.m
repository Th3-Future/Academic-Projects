function [v1,v2,STM,f1,f2,r2sol]=pcrtbplvl1(r1,r2,t1,t2,v1g)

%this is a level 1 differential corrector for a multiple shooting method

%initialize free parameter

e1=[1;0;0;0];
e2=[0;1;0;0];
e3=[0;0;1;0];
e4=[0;0;0;1];

zfree=v1g;

zfix=r1;
T2i=[zeros(2,2) eye(2);eye(2) zeros(2,2)];
T2ifree=[zeros(2,2);eye(2)];
T1=[eye(2) zeros(2,2)];
Tspan=[t1 t2];

pdes=r2;

options=odeset('RelTol',1e-13,'AbsTol',1e-17);

for jj=1:10,
    
    z=[zfree;zfix];
    X0=[T2i*z;e1;e2;e3;e4];
    
    [Tsp,Xsp] = ode45('pcrtbp_dyn_stm',Tspan,X0,options);

    %extract final state and state transition matrix at end point

    [ros,~]=size(Tsp);

    Xf=Xsp(ros,1:4).';

    STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];
    
    HH=T1*STM*T2ifree;
    
    zfree=zfree+HH\(pdes-T1*Xf);
    
end
    

z=[zfree;zfix];
X0=[T2i*z;e1;e2;e3;e4];
    
[Tsp,Xsp] = ode45('pcrtbp_dyn_stm',Tspan,X0,options);

%extract final state and state transition matrix at end point

[ros,~]=size(Tsp);

Xf=Xsp(ros,1:4).';

v1=zfree;

v2=Xf(3:4,1);

STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];
f1=pcrtbp_dyn(t1,X0(1:4,1));
f2=pcrtbp_dyn(t1,Xf);

r2sol=Xf(1:2,1);


end
