from math import*
import Cl_fetch
import time
start=time.time()
teta_list=[18.53,12.79,9.25,6.87,5.16,3.88,2.89,2.09,1.44,0.9,0.44];
r_list=[0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3] #design variables to be obtained
c_list=[0.482,0.369,0.297,0.248,0.213,0.186,0.166,0.149,0.135,0.124,0.114];F=0.85
n=0
ilist=[['r[m]','teta[deg]','tsr_r','c[m]','phi[deg]','alpha[deg]','Cl','a',
        'an']]
while n<len(r_list):
    def zero(f,lb,ub):
        x=lb
        fn=f(x)
        mul=f(lb)*f(ub)
        #print(mul)
        if mul<0:
            x_l=[[],[]]
            i=0
            while round(fn,4)!=0 and x<ub:
                fn=f(x)
                x+=0.0001
                x_l[0].append(x)
                x_l[1].append(fn)
                #print(x_l[0][i], x_l[1][i])
                i+=1
                if fn>0.0001 and round(fn,4)!=0:
                    break
            if not(round(fn,4)==0):
                n=min([n for n in x_l[1] if n>0])
                j=x_l[1].index(n)
                x=x_l[0][j]
                fn=x_l[1][j]
                print('when r is very small, phi',x )
            return x,fn  
        else:
            null='null'
            return null
        
    teta=radians(teta_list[n]); r=r_list[n]; c=c_list[n]
    B=3;  r_tip=3; tsr=7; Vd=7; V=7

    omega = round((tsr * Vd) / r_tip, 4)
    tsr_r = round((omega * r) / V, 4)
    phi = (2/3)*(atan(1/tsr_r))
    #c = round((8*pi*r*sin(phi))/(3*B*tsr_r), 4)
    alpha = round(degrees(phi - teta), 2)
    def Cl_search(alpha):
        return round(Cl_fetch.Cl_search(alpha), 4)
    Cl=0.85#Cl_search(alpha)
    print('\nradius',r,' phi',round(degrees(phi),2), 'alpha',alpha)
    sigma = round((B * c)/(2* pi * r), 5)
    print('Chord lenght',c,'Cl',Cl,'sigma=',sigma,' tsr_r',tsr_r,'omega',omega)
    def k(phi,Cl,F):
        return ((sigma*Cl*cos(phi))/(4*F*sin(phi)**2))
    def k_t(phi,Cl,F):
        return ((sigma*Cl)/(4*F*cos(phi)))
    def kv(k):
        return (k/(1+k))
    def y(k,F):
        y1=(2*F*k-(10/9-F)); y2=(2*F*k-F*(4/3-F)); y3=(2*F*k-(25/9-2*F))
        return (y1-sqrt(y2))/y3
    def axial(k,y):
        if k>-1 and k<2/3:
            chk=1
            print('Momentum',k/(1+k))
            return (k/(1+k)), chk
        if k>=2/3:
            chk=2
            print('Emperical',y)
            return y, chk
        else:
            print('No solution')
            return 'no_soln'
    def angular(k_t):
        return (k_t/(1-k_t)) 
    rst=axial(k(phi,Cl,F), y(k(phi,Cl,F),F))
    a=rst[0]; chk=rst[1]
    if a=='no_soln':
        break
    #an=((1 - 3*a)/(4*a-1))
    an=angular(k_t(phi,Cl,F))
    #print('an',an,'k_t',k_t(phi,Cl,F),'phi',degrees(phi))
    e=1*10**-6
    def fn(phi):
        return ((sin(phi)/(1-a))-(cos(phi)/tsr_r*(1+an)))
    def fpb(phi):
        kn=k(phi,Cl,F)
        return ((sin(phi)*(1-kn))-(cos(phi)/tsr_r*(1+an)))
    if fn(pi/2)>0 and not(an<-1):
        soln_new = zero(fn,e,pi/2)
        phi_new = soln_new[0]
        #print('inflow angle',round(degrees(phi_new), 2),'[deg]')
    elif fpb(-pi/4)<0 and fpb(e)>0:
        print('Propeller brake region')
        soln_new = zero(fpb,-pi/4,-e)
        phi_new = soln_new[0]
    else:
        print('locally reversed tangential flow')
        soln_new = zero(fn,pi/2,pi)
        phi_new = soln_new[0]
    alpha = round(degrees(phi_new - teta), 2)
    Cl = Cl_search(alpha)
    if chk==1:
        a_f = kv(k(phi_new,Cl,F))
    else:
        a_f = y(k(phi_new,Cl,F),F)
    #an_f = angular(k_t(phi_new,Cl,F))
    an_f=((sigma*Cl)/(4*tsr_r*sin(phi_new))*(1-a_f))
    def fn2():
        return (F*(tsr_r**3)*an_f*(1-a_f))
    fn_tsr_r=fn2()
    ilist.append([r,degrees(teta),tsr_r,c,degrees(phi_new),alpha,Cl,a_f,an_f,
                  fn_tsr_r])
    n+=1


dash='-'*80
for i in range(len(ilist)): #tabular printing
    if i==0:
        print(dash)
        print('{:>0s}''{:>10s}''{:>8s}''{:>9s}''{:>12s}''{:>11s}''{:>5s}'
              '{:>8s}''{:>10s}'
              .format(ilist[i][0],ilist[i][1],ilist[i][2],ilist[i][3],ilist[i][4],
                      ilist[i][5],ilist[i][6],ilist[i][7],ilist[i][8]))
        print(dash)
    else:
        print(('{:>0.1f}''{:>9.2f}''{:>10.3f}''{:>10.3f}''{:>9.2f}''{:>9.2f}'
               '{:>10.4f}''{:>10.4f}''{:>10.4f}'
              .format(ilist[i][0],ilist[i][1],ilist[i][2],ilist[i][3],ilist[i][4]
                      ,ilist[i][5],ilist[i][6],ilist[i][7],ilist[i][8])))
fn_tsrr_s=0
for i in range(len(r_list)-1):
    tsrr_m=((ilist[i+2][2]-ilist[i+1][2])/2)
    fn_exp=((ilist[i+1][9]+ilist[i+2][9]))
    fn_tsrr_m=(tsrr_m)*fn_exp
    fn_tsrr_s+=fn_tsrr_m
    #print('mean tsr_r',tsrr_m,'function expand',fn_exp,'each integra f(tsr_r)', fn_tsrr_m)
Cp=(8/tsr**2)*(fn_tsrr_s)
##for i in range(0,len(ilist)-1):
##    print('fn=',ilist[i+1][9])
#print('k', k(phi,Cl,F), 'k_new', k(phi_new,Cl,F))
print('\nsum integra f(tsr_r)',fn_tsrr_s, 'Cp',Cp,)
end=time.time()
time=end-start 
print('\nsimulation time:',round(time),'[sec]')


