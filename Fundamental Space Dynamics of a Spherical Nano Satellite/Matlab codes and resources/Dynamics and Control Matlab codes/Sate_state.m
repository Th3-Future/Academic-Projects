function dstatedt = Sate_state(t,state)
global BB invI I m nextMagUpdate lastMagUpdate lastSensorUpdate 
global nextSensorUpdate BfieldMeasured pqrMeasured BfieldNav pqrNav
global BfieldNavPrev pqrNavPrev
x = state(1);
y = state(2);
z = state(3);

q0123 = state(7:10);
p = state(11);
q = state(12);
r = state(13);
pqr = state(11:13);
%% Kinematics
vel = state(4:6);% v = [xdot;ydot;zdot]
%% Rotational Kinematics
PQRMAT = [0 -p -q -r;p 0 r -q;q -r 0 p;r q -p 0];
q0123dot = 0.5*PQRMAT*q0123;
%% Earth Planet Parameters
mu = 3.986e14; %m3/s2
%% Gravity Model
r_e = state(1:3); % r = [x;y;z]
rho = norm(r_e);
rhat = r_e/rho;
Fgrav = -(mu*m/rho^2)*rhat;
%%% Dynamics
accel = Fgrav/m;
%% The magnetic field model
if t >= lastMagUpdate
    lastMagUpdate = lastMagUpdate + nextMagUpdate;
    %%%Convert Cartesian x,y,z into Lat,Lon, Alt
    phiE = 0;
    thetaE = acos(z/rho);
    psiE = atan2(y,x);
    latitude = 90-thetaE*180/pi;
    longitude = psiE*180/pi;
    rhokm = (rho)/1000;
    [BN,BE,BD] = igrf('01-Jan-2020',latitude,longitude,rhokm,'geocentric');
    %%%Convert NED (North East Down to X,Y,Z in ECI frame)
    %%%First we need to create a rotation matrix from the NED frame to the 
    %%%inertial frame
    BNED = [BN;BE;BD]; 
    BI = TIB(phiE,thetaE+pi,psiE)*BNED;
    %BI = eye(3)*BNED;    
    BB = TIBquat(q0123)'*BI;
    %%%Convert to Tesla
    BB = BB*1e-9;
end
if t >= lastSensorUpdate
    %%%%SENSOR BLOCK
    lastSensorUpdate = lastSensorUpdate + nextSensorUpdate;
    [BfieldMeasured,pqrMeasured] = Sensor(BB,pqr); 
    
    %%%NAVIGATION BLOCK
    [BfieldNav,pqrNav] = Navigation(BfieldMeasured,pqrMeasured);   
end


%% CONTROL BLOCK
current = Control(BfieldNav,pqrNav);
%%%magtorquer parameters
n = 84;
A = 0.02;
muB = current*n*A;

%%%Magtorquer Model
LMN_magtorquers = cross(muB,BB);
%% Rotational Dynamics
C=0.00003;
H = I*pqr;
%Cf=(rp^3*pqr.^2)/mu; %coefficient of friction
pqrdot = invI*((LMN_magtorquers-pqr.*C) - cross(pqr,H));

%% Return derivatives vector
dstatedt = [vel;accel;q0123dot;pqrdot];

