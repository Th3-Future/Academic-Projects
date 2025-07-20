!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves an unsteady linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS1 q4a
!****************************End Header*****************************************
program unsteady_Vfinite
    implicit none
    !Variable declaration
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: Ti, b_nth
    DOUBLE PRECISION, PARAMETER :: lamda=12, rho=8933, Cp=385, h=50, L=0.25
    DOUBLE PRECISION :: T_wall, dt, tt
    INTEGER :: N
    N=10*25 !Number of Control Volume
    dt=16.91     !Our choosen timestep (del_t)
    ALLOCATE(Ti(n),b_nth(n))
    Ti=273.15        !Initial Temperature at t=0
    T_wall=273.15     !Initial Temperature at x=0 and t=0
    call unsteady(b_nth,Ti,T_wall,tt) !solves the unsteady system
    
    Print *, "t_total = ", (tt-dt)*1D-5,"sec"
    
    contains
!*************************************************************************************
    subroutine unsteady(b, Tp, Tww, t)
        DOUBLE PRECISION, dimension(n) :: b, e, f, g, Tp, xi
        DOUBLE PRECISION :: aa, dx, r, Tww, alpha, t
        integer  :: i, j
    !!This subroutine solves unsteady thermal problem using O'Brien's method

        open(1, file = 'PS2_Q4a.txt', status = 'unknown') 
        !Variable initialization
        t=0; e=0; g=0
        dx=L/n
        r=dt/dx**2
        alpha=lamda/(rho*Cp) !Thermal diffusivity
        aa=h*dx/(2*lamda)
        do while (Tww<0.99*298.15) !Loop for time step
            Tww=Tp(1)+aa*(298.15-Tww) !Temperature at the surface from our formulation
            do i = 1,n !Loop for control volume
                if ( i==1 ) then
                    g(1)=-alpha*r
                    f(1)=1+alpha*r
                    b(1)=alpha*dt*h*(298.15-Tww)/(lamda*dx)+Tp(1) !The RHS of our formulation
                    xi(1)=dx/2
                else if ( i<n ) then
                    g(i)=-alpha*r
                    f(i)=1+2*alpha*r
                    e(i)=-alpha*r
                    b(i)=Tp(i)
                    xi(i)=xi(i-1)+dx
                else
                    f(n)=1+alpha*r
                    e(n)=-alpha*r
                    b(n)=Tp(n)
                    xi(i)=xi(i-1)+dx
                end if  
            end do
            !For computational efficiency, the A tri-diagonal matrix is split to 3 vectors
            call tdma(e,f,g,b,Tp) !Solves the linear system at each time-step
            if ( t==0 ) then
                write(1,3) "Time(sec)", "T(0,t)", (xi(j), j=1,n) !writes the result's title to .txt file
            end if
            write(1,4) t*1D-5, Tww, (Tp(j), j=1,n) !writes the results to .txt file
            t=t+dt
        end do
        3 format(a9,3x,a6,250(2x,f6.3))
        4 format(f9.7,3x,f7.3,250(1x,f7.3))
        close(1)
    end subroutine unsteady
!*************************************************************************************
!*************************************************************************************
    subroutine tdma(e, f, g, b1, x)
        DOUBLE PRECISION, dimension(n), INTENT(IN):: b1
        DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: x
        DOUBLE PRECISION, DIMENSION(n):: b, e, f, g
        INTEGER :: k
    !!This subroutine solves a tri-diagonal linear system using the Thomas Algorithm
    
        b=b1
       !Decomposition
        do k = 2,n
            e(k) = e(k)/f(k-1)
            f(k) = f(k) - e(k)*g(k-1)
        end do
       !Forward Substitution
        do k = 2,n
            b(k) = b(k) - e(k)*b(k-1)
        end do
       !Backward Substitution
        x(n)=b(n)/f(n)
        do k = n-1,1,-1  !(step size of –1)
            x(k) = (b(k) - g(k)*x(k+1))/f(k)
        end do
    end subroutine tdma
!*************************************************************************************
end program unsteady_Vfinite