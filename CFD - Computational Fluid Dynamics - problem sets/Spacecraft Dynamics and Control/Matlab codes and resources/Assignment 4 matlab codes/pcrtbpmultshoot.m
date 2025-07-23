function [X0,t0,tf,PP,VV,YY]=pcrtbpmultshoot(PP,Vg)

%this function implements a multiple shooting method for the PCRTBP
%dynamics

%PP is an array containing the position vectors and epochs of all of n patch points,
%with each column having the form [R;t], where R is position, and t the
%epoch

%Vg is an array containing an initial guess at the velocity vectors for the
%first n-1 patch points

[~,N]=size(PP); %number of patch points

Vi=Vg; %initialize V1+ for each segment

STMs=zeros(4,4,N-1); %initialize 3 dimensional array of state transition matrices for each of the segments

VV=zeros(2,N-1,2); %initialize 3 dimensional array of end point velocities for each segment

ff=zeros(4,N-1,2); %initialize 3 dimensional array of end point dynamics for each segment



for jj=1:10,

%Level 1 corrector (determine V1+ and V2- for each segment of the
%trajectory)

    for kk=1:N-1,
        [v1,v2,STM,f1,f2,~]=pcrtbplvl1(PP(1:2,kk),PP(1:2,kk+1),PP(3,kk),PP(3,kk+1),Vi(:,kk));
        VV(:,kk,1)=v1;
        VV(:,kk,2)=v2;
        STMs(:,:,kk)=STM;
        ff(:,kk,1)=f1;
        ff(:,kk,2)=f2;   
    end

    
    %Level 2 corrector
    %form column of all Delta Vs.
    YY=zeros(2*(N-2)+5,1);  %initialize vector containing all Delta Vs and additional constraints
    pp=0;
    for kk=2:N-1,
        YY((2*pp+1):(2*pp+2),1)=VV(:,kk,1)-VV(:,kk-1,2);
        pp=pp+1;
    end
    
    M=zeros(2*(N-2)+5,3*N);   %initialize Jacobian
    %form Jacobian
    pp=0;
    for kk=1:N-2,
        STM1=STMs(:,:,kk);
        STM2=STMs(:,:,kk+1);
        
        STM1rr=STM1(1:2,1:2);
        STM1rv=STM1(1:2,3:4);
        STM1vr=STM1(3:4,1:2);
        STM1vv=STM1(3:4,3:4);
        
        STM2rr=STM2(1:2,1:2);
        STM2rv=STM2(1:2,3:4);
        STM2vr=STM2(3:4,1:2);
        STM2vv=STM2(3:4,3:4);
        
        STM1r=[STM1rr STM1rv];
        STM1v=[STM1vr STM1vv];
        
        STM2r=[STM2rr STM2rv];
        
        f11=ff(:,kk,1);
        f12=ff(:,kk,2);
        
        f21=ff(:,kk+1,1);
        f22=ff(:,kk+1,2);
        
        M1=-STM1vr+STM1vv*(STM1rv\STM1rr);
        M2=STM1v*f11-STM1vv*(STM1rv\(STM1r*f11));
        M3=-STM2rv\STM2rr-STM1vv/STM1rv;
        M4=STM2rv\(STM2r*f21)-f12(3:4,1)+STM1vv*(STM1rv\f12(1:2,1));
        M5=STM2rv\eye(2);
        M6=-STM2rv\f22(1:2,1);
        
        M((2*pp+1):(2*pp+2),(3*pp+1):(3*pp+9))=[M1 M2 M3 M4 M5 M6];
        
        pp=pp+1;
    end
    
    %additional output and Jacobian rows for constraints
    %initial time constraint
    YY(2*(N-2)+1,1)=PP(3,1);
    M(2*(N-2)+1,3)=1;
    
    %initial y-coordinate constraint
    YY(2*(N-2)+2,1)=PP(2,1);
    M(2*(N-2)+2,2)=1;
    
    %initial x-velocity constraint
    YY(2*(N-2)+3,1)=VV(1,1,1);
    
    STM1=STMs(:,:,1);
    
    STM1rr=STM1(1:2,1:2);
    STM1rv=STM1(1:2,3:4);

    STM1r=[STM1rr STM1rv];

    f11=ff(:,1,1);
    f12=ff(:,1,2);
    M1=[1 0]*(-STM1rv\STM1rr);
    M2=[1 0]*(STM1rv\(STM1r*f11));
    M3=[1 0]*(STM1rv\eye(2));
    M4=[1 0]*(-STM1rv\f12(1:2,1));
    M(2*(N-2)+3,1:6)=[M1 M2 M3 M4];
    
%     %final y-position constraint
%     YY(2*(N-2)+4,1)=PP(2,N);
%     M(2*(N-2)+4,3*N-1)=1;
%     
%     %final x-velocity constraint
%     YY(2*(N-2)+5,1)=VV(1,N-1,2);
%     
%     STM1=STMs(:,:,N-1);
%     
%     STM1rr=STM1(1:2,1:2);
%     STM1rv=STM1(1:2,3:4);
%     STM1vr=STM1(3:4,1:2);
%     STM1vv=STM1(3:4,3:4);
% 
%     STM1r=[STM1rr STM1rv];
%     STM1v=[STM1vr STM1vv];
% 
%     f11=ff(:,N-1,1);
%     f12=ff(:,N-1,2);
%     M1=[1 0]*(STM1vr-STM1vv*(STM1rv\STM1rr));
%     M2=[1 0]*(-STM1v*f11+STM1vv*(STM1rv\STM1r*f11));
%     M3=[1 0]*(STM1vv/STM1rv);
%     M4=[1 0]*(f12(3:4,1)-STM1vv*(STM1rv\f12(1:2,1)));
%     M(2*(N-2)+5,(3*N-5):3*N)=[M1 M2 M3 M4];

    %final position constraint
    YY((2*(N-2)+4):(2*(N-2)+5),1)=PP(1:2,1)-PP(1:2,N);
    M((2*(N-2)+4):(2*(N-2)+5),1:2)=eye(2);
    M((2*(N-2)+4):(2*(N-2)+5),(3*N-2):(3*N-1))=-eye(2);
    
%     %final velocity constraint
%     YY((2*(N-2)+6):(2*(N-2)+7),1)=VV(1:2,1,1)-VV(1:2,N-1,2);
%     STM1=STMs(:,:,1);
%     
%     STM1rr=STM1(1:2,1:2);
%     STM1rv=STM1(1:2,3:4);
% 
%     STM1r=[STM1rr STM1rv];
% 
%     f11=ff(:,1,1);
%     f12=ff(:,1,2);
%     M1=(-STM1rv\STM1rr);
%     M2=(STM1rv\(STM1r*f11));
%     M3=(STM1rv\eye(2));
%     M4=(-STM1rv\f12(1:2,1));
%     M((2*(N-2)+6):(2*(N-2)+7),1:6)=[M1 M2 M3 M4];
%     
%     STM1=STMs(:,:,N-1);
%     
%     STM1rr=STM1(1:2,1:2);
%     STM1rv=STM1(1:2,3:4);
%     STM1vr=STM1(3:4,1:2);
%     STM1vv=STM1(3:4,3:4);
% 
%     STM1r=[STM1rr STM1rv];
%     STM1v=[STM1vr STM1vv];
% 
%     f11=ff(:,N-1,1);
%     f12=ff(:,N-1,2);
%     M1=(STM1vr-STM1vv*(STM1rv\STM1rr));
%     M2=(-STM1v*f11+STM1vv*(STM1rv\STM1r*f11));
%     M3=(STM1vv/STM1rv);
%     M4=(f12(3:4,1)-STM1vv*(STM1rv\f12(1:2,1)));
%     M((2*(N-2)+6):(2*(N-2)+7),(3*N-5):3*N)=-[M1 M2 M3 M4];
    %patch point correction
    DPPvec=-M.'*((M*M.')\YY);
    DPP=zeros(3,N);
    for kk=0:N-1,
        DPP(:,kk+1)=DPPvec((3*kk+1):(3*kk+3),1);
    end
    YY
    PP=PP+DPP;
    
    %set initial segment velocity guesses for next iteration
    Vi=VV(:,:,1);
    

end

%get patch point input and output velocities for each patch point

for kk=1:N-1,
        [v1,v2,STM,f1,f2,~]=pcrtbplvl1(PP(1:2,kk),PP(1:2,kk+1),PP(3,kk),PP(3,kk+1),Vi(:,kk));
        VV(:,kk,1)=v1;
        VV(:,kk,2)=v2;
end

X0=[PP(1:2,1);VV(:,1,1)];
t0=PP(3,1);
tf=PP(3,N);