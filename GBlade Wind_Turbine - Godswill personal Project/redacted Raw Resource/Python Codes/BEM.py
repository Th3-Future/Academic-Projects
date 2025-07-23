from math import*
import BEM_back
#! /usr/bin/octave -qf
# Octave script for bem calculation
# Grant Ingram (c) 2007

# this does just one radial element - mainly for demonstration purposes...

# blade details
r = 5; gamma = 92.6; r_tip = 5
print("\nElement: r:"+str(r)+"[m]  gamma:" +str(gamma)+"[deg] \n")

# wind details
V=7; B=3; tsr=8
omega = round((tsr * V) / r_tip, 4) # note blade speed is derived
tsr_r = round((omega * r) / V, 4)
beta = round((90 - ((2/3)*degrees(atan(1/tsr_r)))), 2)
c =0.19 #round((8*pi*r*cos(radians(beta)))/(3*B*tsr_r), 4) # chord distribution
sigma = round((B * c) / (2* pi * r), 5)
print("V:"+str(V)+ "[m/s]  Omega:"+str(omega)+ "[rad/s]  tsr_r:"+str(tsr_r)+ ' Chord:'+str(c)+
      '[m]'+'\n  Solidity Sigma:',sigma,'  Intial Beta:',beta,'\n')

# guess initial value
i = round(gamma - beta, 2)
print('i', i)
CL = round(BEM_back.Cl_search(i), 4)
k=0
axial = 1/(1 + ((4 * cos(radians(beta))**2)/(sigma *CL* sin(radians(beta)))))
angular = ((1 - 3*axial)/(4*axial -1))
print("Initial values: axial:", round(axial, 4), "  angular:", round(angular, 4), "\n")

# iterate
error_axial = 10; error_angular = 10
##old_axial = 5; old_angular = 5;
j=0

while abs(error_axial) > 1*10**-5 and abs(error_angular) > 1*10**-5:
    j=1+j
    k=(sigma *CL* sin(radians(beta)))/(4*cos(radians(beta))**2)
    #if k<2/3:
    old_axial = axial; old_angular = angular
    beta = degrees(atan ( (tsr_r * ( 1 + angular)) / ( 1- axial)))
    i = round(gamma - beta, 2) # converted to degrees here
    CL = round(BEM_back.Cl_search(i), 4) # find new values of CL,CD from the aerofoil curve
    axial = (( sigma *( CL*sin(radians(beta)))) / (4*cos(radians(beta))**2))*(1-axial)
    #axial = 1 / ((1 + (4 * cos(radians(beta))**2)/(sigma*CL*sin(radians(beta)))))
    angular = (( sigma * CL ) / (4* tsr_r * cos(radians(beta))))*(1-axial)
    error_axial = axial - old_axial
    error_angular = angular-old_angular
##    else:
##        break
    if j<4:
        print('k',k, 'Beta=',beta)
        print('axial',round(axial,4),'angular',round(angular,4),'Cl',CL)
print('number of itertion',j)
V_r = round((omega*r)/tsr,4)
Vr = round(sqrt((V_r**2)+(omega**2)),4)
c_opt = 'not a local radius'
if r==r_tip:
    c_opt = round((2*pi*r*8*V)/(B*9*CL*tsr_r*Vr),4)
print("Final values: axial:", round(axial,4), "  angular:", round(angular,4), "\n")
print("Flow parameters: i:"+str(i)+ "[deg]  CL:"+str(CL)+ 
      "  beta:"+str(round(beta, 2))+'[deg]')
print('\nlocal wind speed:',V_r,'Relative wind speed:',Vr,'Optimum Chord:',c_opt)
