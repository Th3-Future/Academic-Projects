clear; clc
%% Variable Initialization
q0=[50000,75000];
Tin=18+273.15;   %Temperature at the inlet
t=0.00006;       %wall thickness
L2=0.0374;       %Lenght of heated part 2
L1=0.0272;       %Lenght of part 1
k_s=237;         %Thermal conductivity of aluminium
n=L2/50;         %number of grid point along z-direction
z=(0:n:L2);
zn=z-L2/2;
size=size(z);
col=size(2);
fig1 = figure(1); ax1 = axes; hold(ax1,'on');
%% Loop through the two independent heat flux
for j=1:2
    if j==1
        for i=1:col   %Loops through each grid point for q0_1
            q=q0(j)*(1-(4*zn(i)^2/L2^2)); % from equation (1.1)
            T(i)=(q*t/k_s)+Tin;           % from equation (3.1)
        end
    plot(ax1,(z+L1)*1000,T-273.15)
    else
        for i=1:col   %Loops through each grid point for q0_2
            q=q0(j)*(1-(4*zn(i)^2/L2^2)); % from equation (1.1)
            T(i)=(q*t/k_s)+Tin;           % from equation (3.1)
        end
       
    plot(ax1,(z+L1)*1000,T-273.15)
    %Plot settings
    title('Analytical Solution for Simple Heat Conduction');
    xlabel('Z $(mm)$','Interpreter','latex','fontsize',14);
    ylabel('T $(degC)$','Interpreter','latex','fontsize',14);
    leg1 = legend('$q_0=50000$','$q_0=75000$');
    set(leg1,'Interpreter','latex');
    end
end
