% -------------------------------------------------------------------
%  Written by Godswill Ezeorah for AE8137 Course Project on 
%  13-Nov-2021, to automate Simulink Simulation 
% -------------------------------------------------------------------
clear; clc;
%% Variable Initialization [table (2)]
Ix = 1.2; Iy = 0.5; Iz = 0.9;              % Refrence axis Moment of inertial 
Iry = 0.07;                                % Rotor y-axis Moment of inertial
theta_dot = 9;                             % Rotor angular velocity for raduis=3m)
Kmoi = (Ix-Iy)/Ix*Iy*Iz;
R_tol = 1e-12;                             % Solver Reletive Tolerance
%% Desired Inputs for Tracking [table (1)]
wind = [30;30;30;50;70;90;70;50;30;30];    % Angular position
wind_dot = [0;0;0;20;20;20;-20;-20;-20;0]; % Angular velocity
wind_ddot = [0;0;0;20;0;0;-40;0;0;20];     % Angular accleretion
time_step = [1;2;3;4;5;6;7;8;9;10];        % corresponding time

%% Control Performance Specification
Tsmax = 1; Tsmin = 0.012;               
OS = [0,0,8];                              % range Percentage Overshot 
[r,n]=size(OS);
mT = (Tsmax-Tsmin)/n;
T_s = (Tsmin:mT:Tsmax);                   % range Settling Time
[row,col] = size(OS);

%% Looping Through all Performance Specification in Simulink
for i=1:col
    % Feedback gain definition
    K1 = (16/T_s(i)^2)*(pi^2/log(OS(i)/100)^2+1);     % Equation (2.6)
    K2 = 8/T_s(i) ;                                   % Equation (2.7)
    
    %% Extraction of Simulink Output Dataset
    sim1=sim('yaw_control.slx');               %Gets the simulink simulation dataset
    sig = sim1.yout.getElement('output');  %Extracts time structured output signal
    t = sig.Values.time;                   %Extracts time array
    phi_a = sig.Values.Data(:,2);          %Extracts measured position
    phi_d = sig.Values.Data(:,1);          %Extracts desired position
    torque = sig.Values.Data(:,3);         %Extracts the new torque signal
    
    %% Plotting the Simulation Results
    fig1 = figure(i); ax1 = axes; hold(ax1,'on'); % for Control plot
    plot(ax1,t, phi_d);
    plot(ax1,t, phi_a,'--r');
    leg1 = legend('$phi$','Desired Angle');
    set(leg1,'Interpreter','latex');
    title(['Control System Performance (%OS=',num2str(OS(i)),...
        ' and T_s=',num2str(T_s(i)),')']);
    xlabel('Time $(sec)$','Interpreter','latex','fontsize',12);
    ylabel('Yaw Angle $(deg)$','Interpreter','latex','fontsize',12);
    axis([0 10 0 95]);
    set(gca,'Color',[0.9,0.9,0.9]);
    
    fig2 = figure(i+3); ax2 = axes; hold(ax2,'on'); % for Torque plot
    plot(ax2,t, torque);
        title(['Yaw Control Torque (%OS=',num2str(OS(i)),...
        ' and T_s=',num2str(T_s(i)),')']);
    xlabel('Time $(sec)$','Interpreter','latex','fontsize',12);
    ylabel('Torque $(N.m)$','Interpreter','latex','fontsize',12);
    axis([0 10 -5e3 4.5e3]);
    set(gca,'Color',[0.9,0.9,0.9]);
end
