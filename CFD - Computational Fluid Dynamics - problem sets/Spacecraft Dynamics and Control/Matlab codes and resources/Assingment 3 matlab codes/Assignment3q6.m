clear
clc

mu=0.01213;         %mass ratio for moon

coeff=[-1 mu+2 -2*mu-1 -mu+1 2*mu-2 1-mu];       %coefficients for polynomial to solve for L2 location

C=roots(coeff)

D=imag(C)

D=abs(D)

[~,I]=min(D)

rs1=C(I)

rs2=1-rs1

xL1=rs1-mu

yL1=0


%Linearization about L2 (in-plane only)

A=(1-mu)/rs1^3+mu/rs2^3

F=[0 0 1 0;
    0 0 0 1;
    1+2*A 0 0 2;
    0 1-A -2 0]

%Full state vector at L2
XL2=[xL1;
    yL1;
    0;
    0];

%eigendecomposition of F

[V,D]=eig(F)



%center manifold frequency
om=imag(D(3,3));

%initial guess of half period pi/om
Tp=pi/om

%real and imaginary parts of corresponding eigenvector
uR=real(V(:,3));
uI=imag(V(:,3));

%x-axis crossing

x0=xL1+2*0.001*uR(1,1);

%initial y-axis velocity
dy0=2*0.001*uR(4,1)
options=odeset('RelTol',1e-12,'AbsTol',1e-15);
T1=[0 1 0 0;
    0 0 1 0;
    1 0 0 0;
    0 0 0 1];
T2=[0 0 0 1;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0];
%set initial conditions (stack columns of state transition matrix under
%state vector)
X0=[x0;0;0;dy0;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];

for jj=1:7,

[Tsp,Xsp] = ode45('pcrtbp_dyn_stm',[0 Tp],X0,options);

%extract final state and state transition matrix at end point

[ros,~]=size(Tsp);

Xf=Xsp(ros,1:4).';

STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];

%correction
Am=T1*STM*T2.';
Amp=Am(1:2,1);
dX=pcrtbp_dyn_stm(1,Xsp(ros,:).');
Em=T1*dX(1:4,1);
Emp=Em(1:2,1);

Q=[Amp Emp];

dpar=Q\[-Xf(2,1);-Xf(3,1)];

ddy0=dpar(1,1);

dTp=dpar(2,1);

dy0=dy0+ddy0;

Tp=Tp+dTp;

X0=[x0;0;0;dy0;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];

end
% [Tsp,Xsp] = ode45('pcrtbp_dyn_stm',[0 2*Tp],X0,options);
% 
% %extract final state and state transition matrix at end point
% 
% [ros,~]=size(Tsp);
% 
% Xf=Xsp(ros,1:4).';
% 
% STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];
% 
% eig(STM)
%converged period
Tp2=Tp*2;

%converged initial condition
X0=[x0;0;0;dy0];

%final periodic orbit

[Tsp,Xsp] = ode45('pcrtbp_dyn',[0 Tp2],X0,options);


figure(1)
%plot(Xsp(:,1),Xsp(:,2))
hold on
plot(1-mu,0,'ko')
plot(xL1,0,'r+')
title('Lyapunov Orbit Family About L_1')
leg1=legend('Moon','L_1')



%compute eigenvalues of the monodromy matrix for z-direction
X0=[x0;0;0;dy0;1;0;0;1];
    [Tsp,Xsp] = ode45('pcrtbp_dyn_stmz',[0 Tp2],X0,options);
    %extract state transition matrix at end point

    [ros,~]=size(Tsp);

    STM=[Xsp(ros,5:6).' Xsp(ros,7:8).'];

    Dstm=eig(STM)

figure(2)
hold on
plot(real(Dstm),imag(Dstm),'k.');
title('Location of Bifurcation into new Periodic Orbit Families from L_1')

figure(3)
hold on
plot(0.001,max(abs(Dstm)),'k.')
title('Direction in which these new Families Grows')
kk=0;
%expand family of periodic orbits, with fixed x-axis crossings.
for rho=0.002:0.001:0.002+120*0.001,
    x0=xL1+2*rho*uR(1,1);
    X0=[x0;0;0;dy0;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
    kk=kk+1;
    for jj=1:7,

    [Tsp,Xsp] = ode45('pcrtbp_dyn_stm',[0 Tp],X0,options);

    %extract final state and state transition matrix at end point

    [ros,~]=size(Tsp);

    Xf=Xsp(ros,1:4).';

    STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];

    %correction
    Am=T1*STM*T2.';
    Amp=Am(1:2,1);
    dX=pcrtbp_dyn_stm(1,Xsp(ros,:).');
    Em=T1*dX(1:4,1);
    Emp=Em(1:2,1);

    Q=[Amp Emp];

    dpar=Q\[-Xf(2,1);-Xf(3,1)];

    ddy0=dpar(1,1);

    dTp=dpar(2,1);

    dy0=dy0+ddy0;

    Tp=Tp+dTp;

    X0=[x0;0;0;dy0;1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];

    end 
    
    if kk==1;
    %converged period
    Tp2=Tp*2

%     %converged initial condition
%     X0=[x0;0;0;dy0];

    %final periodic orbit

    [Tsp,Xsp] = ode45('pcrtbp_dyn_stm',[0 Tp2],X0,options);

    figure(1)
    plot(Xsp(:,1),Xsp(:,2))
    set(leg1,'visible','off')
    
%     %compute eigenvalues and eigenvectors of the monodromy matrix for planar dynamics
%     [ros,~]=size(Tsp);
% 
%     Xf=Xsp(ros,1:4).';
% 
%     STM=[Xsp(ros,5:8).' Xsp(ros,9:12).' Xsp(ros,13:16).' Xsp(ros,17:20).'];
%     
%     dX=pcrtbp_dyn_stm(1,X0);
%     dX(1:4,1)
%     [V,D]=eig(STM)
%     D=(STM-eye(4))*(STM-eye(4));
%     D(1,1)/D(1,4)
%     D(2,1)/D(2,4)
%     D(3,1)/D(3,4)
%     D(4,1)/D(4,4)
    
    %compute eigenvalues of the monodromy matrix for z-direction
    X0=[x0;0;0;dy0;1;0;0;1];
    [Tsp,Xsp] = ode45('pcrtbp_dyn_stmz',[0 Tp2],X0,options);
    %extract state transition matrix at end point

    [ros,~]=size(Tsp);

    STM=[Xsp(ros,5:6).' Xsp(ros,7:8).'];
    
    Dstm=eig(STM);
    figure(2)
    hold on
    plot(real(Dstm),imag(Dstm),'k.');
    figure(3)
    plot(rho,max(abs(Dstm)),'k.')
    kk=0;
    end
end

% %evaluate Jacobi constant along trajectory
% 
% Nsp=length(Tsp);
% 
% Csp=zeros(Nsp,1);
% 
% for ii=1:Nsp,
%     %distances to primaries
%     rs1h=sqrt((Xsp(ii,1)+mu)^2+Xsp(ii,2)^2);
%     rs2h=sqrt((Xsp(ii,1)-1+mu)^2+Xsp(ii,2)^2);
%     Csp(ii,1)=2*((Xsp(ii,1)^2+Xsp(ii,2)^2)/2+(1-mu)/rs1h+mu/rs2h)-Xsp(ii,3)^2-Xsp(ii,4)^2;
% end
% 
% 
% 
% figure(3)
% plot(Tsp,Csp)

%Lyapunov orbit close to bifurcation with Halo orbit family

rho=0.114;

%find 



