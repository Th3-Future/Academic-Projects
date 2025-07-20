!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on July 07, 2021.
!This program solves a Navier-Stokes equation using finite difference method
!and was written as a solution to AE8112 PS4 q2
!****************************End Header*****************************************
program Navier_Stokes
    implicit none
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ri, P, U, Uo
    DOUBLE PRECISION, PARAMETER :: h=0.5, rin=0.01, rout=0.11, Pout=0, Patm=101325
    DOUBLE PRECISION, PARAMETER :: rho=1.2, mu=1.84D-5, Pi=4*atan(1.0)
    DOUBLE PRECISION :: mdot, t, dt, dr, r, rp, rn, aa, ab, ac, ad, u1, uo1
    DOUBLE PRECISION :: rh
    INTEGER :: n, i
    open(1, file = 'PS4_Q2c_withMu.txt', status = 'unknown') 
    !Variable initialization and definition
    t=0.01
    n=100 !Nuumber of control volumes
    ALLOCATE(ri(n), P(n), U(n), Uo(n))
    dr = (rout-rin)/n
    dt = 0.01  !time step
    Uo = 0
    
    do while (t <= 1)
        t = t+dt
        r=rin
        mdot = (10D-5)*(exp(2*t)-1)
        do i = 1, n !This loop solves velocity using our Contuinity formulation
            rp = r+dr
            rh = r+dr/2
            if (i==1) then   !for first grid point
                U(i)=mdot*(1/dr-1/rh)/((2*pi*h*(rin-dr/2))*(1/rh+1/dr))
            else             !for all other grid points
                U(i)=U(i-1)*(1/dr-1/rh)/(1/rh+1/dr)
            end if
            ri(i)=r
            r=r+dr
        end do
        r=rout
        do i = n, 1, -1 !This loop solves pressure backward using our momentum formulation
            rn = r-dr
            rh = r-dr/2
            aa = (2*rho*dr)*((U(i)+U(i-1))/2-(Uo(i)+Uo(i-1))/2)/dt
            ab = rho*((Uo(i)+Uo(i-1))/2)*(Uo(i)-Uo(i-1))
            ac = mu*(Uo(i)/r-Uo(i-1)/rn)
            ad = 2*mu*(Uo(i-1)-(Uo(i)+Uo(i-1))+Uo(i))/dr
            if ( i == n ) then     !for last grid point
                P(i) = aa+ab-ac-ad+Pout
            elseif ( i >= 2) then  !for all other grid points
                P(i) = aa+ab-ac-ad+P(i+1)
            else                   !for first grid point
                u1=(10D-5)*(exp(2*t)-1)*(2*pi*h*(rin-dr/2))
                uo1=(10D-5)*(exp(2*(t-dt))-1)*(2*pi*h*(rin-dr/2))
                aa = (2*rho*dr)*((U(i)+U1)/2-(Uo(i)+Uo1)/2)/dt
                ab = rho*((Uo(i)+Uo1)/2)*(Uo(i)-Uo1)
                ac = mu*(Uo(i)-Uo1)/rh
                ad = 2*mu*(Uo1-(Uo(i)+Uo1)+Uo(i))/dr
                P(i) = aa+ab-ac-ad+P(i+1)
            end if
            r=r-dr
        end do
        Uo = U
    end do
    
    !Printing results
    do i = 1, n
        write(1,*) ri(i), P(i)+Patm
    end do
    print *, 'Velo. Diff',U(1)-U(n)
    print *, 'press. Diff',P(1)-P(n)
    write(*,2) 'Radius =', ri
    write(*,1) '(Pressure - Patm) =', P
    write(*,2) 'Velocity =', U
    1 format(a20,100E10.3)
    2 format(a11,100f7.3)
    close (1)
end program Navier_Stokes