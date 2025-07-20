!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on July 07, 2021.
!This program solves a Navier-Stokes equation using finite volume method
!and was written as a solution to AE8112 PS4 q2
!****************************End Header*****************************************
program Navier_Stokes
    implicit none
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ri, P, U, Uo, Po
    DOUBLE PRECISION, PARAMETER :: h=0.5, rin=0.01, rout=0.11, Pout=0
    DOUBLE PRECISION, PARAMETER :: rho=1.2, mu=1.84D-5, Pi=4*atan(1.0)
    DOUBLE PRECISION :: mdot, t, dt, dr, r, rp, rn, aa, ab, ac, ad, u1, uo1
    DOUBLE PRECISION :: rph, rnh, L2, Tol
    INTEGER :: n, i, j
    open(1, file = 'PS4_Q2a.txt', status = 'unknown') 
    !Variable initialization and definition
    n=8 
    L2=1
    Tol=1D-4
    j=1
    !The below loop increases the number of gridpoint until tolerance is met.
    do while (L2>Tol) 
        t=0.01
        ALLOCATE(ri(n), P(n), U(n), Uo(n), Po(j))
        dr = (rout-rin)/n
        dt = 0.01  !time step
        Uo = 0
        do while (t <= 1)
            t = t+dt
            r=rin
            mdot = (10D-5)*(exp(2*t)-1)
            do i = 1, n !This loop solves velocity using our contuinity formulation
                rp = r+dr
                if (i==1) then !for first grid point
                    U(i)=mdot/(2*pi*h*rp)
                else           !for all other grid points
                    U(i)=U(i-1)*r/rp
                end if
                ri(i)=r
                r=r+dr
            end do
            r=rout
            do i = n, 1, -1 !This loop solves pressure backward using our momentum formulation
                rn = r-dr
                rnh = r-dr/2
                rph = r+dr/2
                aa = (rho*dr)*(U(i)-Uo(i))/dt
                ab = rho*(r*Uo(i)**2-rn*Uo(i-1)**2)/r
                ac = mu*rph*(Uo(i+1)-Uo(i))/(dr*r)
                ad = mu*rnh*(Uo(i)-Uo(i-1))/(dr*r)
                if ( i == n ) then !for last grid point
                    ac = mu*(Uo(i)-Uo(i))/dr
                    P(i) = aa+ab-ac+ad+Pout
                elseif ( i >= 2 ) then !for all other grid points
                    P(i) = aa+ab-ac+ad+P(i+1)
                else                   !for first grid point
                    u1=(10D-5)*(exp(2*t)-1)/(2*pi*h*rin)
                    uo1=(10D-5)*(exp(2*(t-dt))-1)/(2*pi*h*rin)
                    aa = (rho*dr)*(U(i)-Uo(i))/dt
                    ab = rho*(r*Uo(i)**2-rin*Uo1**2)/r
                    ac = mu*rph*(Uo(i+1)-Uo(i))/(dr*r)
                    ad = mu*r*(Uo(i)-Uo1)/(dr*r)
                    P(i) = aa+ab-ac+ad+P(i+1)
                end if
                r=r-dr
            end do
            Uo = U
        end do
        Po(j)=P(1)
        L2=norm2(Po)/j
        print *, 'No. of gridpoints.',n
        print *, 'L2 Error',L2
        write(1,*) n, L2 
        if ( L2>Tol ) then
            DEALLOCATE(ri, P, U, Uo, Po)
            n=n+n
            j=j+1
        end if
    end do
    close (1)
end program Navier_Stokes