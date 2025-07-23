from math import*
B=3;  r_tip=5; tsr=5; V=7; gamma = 92.6; r = 5
teta=90-gamma
omega = round((tsr * V) / r_tip, 4)
tsr_r = round((omega * r) / V, 4)
phi = (2/3)*(atan(1/tsr_r))
beta = (pi/2-phi)
c = 0.19
#alpha = round(degrees(phi - teta), 2)
Cl = 1.4875
sigma = round((B * c) / (2* pi * r), 5)
def k(phi,Cl):
    return ((sigma*Cl*cos(phi))/(4*sin(phi)**2))
####def k_t(phi,Cl):
####    return (sigma*Cl*sin(phi))/(4*sin(phi)*cos(phi))
def ax(k):
    return k/(1+k)
####def angular(k_t):
####    return k_t/(1-k_t)
#if k(phi,Cl)<2/3:
a=ax(k(phi, Cl))
##    #an=angular(k_t(phi,Cl))
print(k(phi, Cl), 2/3)

#a = ((sigma*Cl*cos(phi))/(4*sin(phi)**2))/(1+((sigma*Cl*cos(phi))/(4*sin(phi)**2)))


#i = round(gamma - beta, 2)
axial = 1/(((4*(cos(beta)**2))/(sigma*Cl*sin(beta)))+1)
#axial1 = (( sigma *( Cl*sin(radians(beta)))) / (4*cos(radians(beta))**2))*(1-axial)
#angular = (1 - 3*axial)/(4*axial -1)
print('phi',round(degrees(phi),2),'beta',round(degrees(beta),2))
print('a', a, '\na1', axial)

