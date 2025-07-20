!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS2 q1c
!****************************End Header*****************************************
program newt_Jac_TDMA
    implicit none
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: T_bar, f_T, Ti, ej, fj, gj
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE :: Jac
    DOUBLE PRECISION, PARAMETER :: Pi = 4*atan(1.0), r=0.10
    DOUBLE PRECISION :: tol, norm, Pin, Pout, Pint, HL_error, Teg, start, finish
    INTEGER, dimension(5):: N_cvs
    INTEGER :: n, i1, tim, count
    N_cvs=(/500,1000,3000,8000,10000/) !Array of control Volumes to be computed
    !Initialization of variables
    tol=10.0**(-8) !The given Newton Tolerance
    Teg=800        !Temperature at the edge of the disc

    do i1 = 1, size(N_cvs) !Loop through for different control volumes in the N_cvs array
        call cpu_time(start) 
        n=N_cvs(i1)
        ALLOCATE(T_bar(n),f_T(n),jac(n,n),Ti(n),ej(n),fj(n),gj(n))
        do tim = 1, 10000 !for timing purpose
            count=1
            norm=1 
            !Initial geuss
            T_bar=Teg
            do while (norm>tol)
                call solve_fT(T_bar,F_T)
                call solve_jx(ej,fj,gj,T_bar) !creates the Tri-diagonal vectors
                call tdma(ej,fj,gj,f_T,Ti)
                T_bar=T_bar+Ti
                norm=norm2(Ti)/n
                count=count+1
            end do
            Pint=pi*0.001*(Teg-T_bar(n))*4*1000 !Another approach to calculating Power-in
            call solve_P(T_bar,Pin,Pout)
            HL_error=100*abs(Pin-Pout)/Pin  !% heat loss error
        end do !for timing purpose
        
        print *, "no. of Control Volume", n
        print *, "no. Iteration before Convergence", count
        write(*,3) "% Heat Loss Error = ", HL_error
        call cpu_time(finish)
        !To get the time taken for one tim_loop run, we divide the time diff. with tim.
        write(*,1) "Program Execution Time = ", ((finish-start)/tim)*1D6,"micro_sec"
        if (n<10000) then
            DEALLOCATE(T_bar,f_T,jac,Ti,ej,fj,gj)
        end if 
    end do
 
    1 format(a40,f8.3,a10)
    3 format(a40,f9.7)

    contains
!*************************************************************************************
    subroutine solve_fT(Tr, fT)
        DOUBLE PRECISION, dimension(n), INTENT(OUT) ::  fT
        DOUBLE PRECISION, dimension(n) :: Tr
        DOUBLE PRECISION :: aa, ab, TW, TP, TE, rr, dr
        INTEGER :: i
    !!This subroutine creates and solve the discretized function

        !Initialization of variables
        dr=r/n
        rr=r
        aa=1000/dr    !for equispaced grid
        ab=(11.34*10.0**(-5))/16
        do i = 1,n
            TP=Tr(i)
            if ( i==1 ) then
                TE=Tr(i+1)
                fT(1)=-aa*(TE-TP) + 2*aa*(TP-Teg) - ab*(rr-dr)*(TE+TP)**4 + ab*16*rr*Teg**4
            else if ( i<n ) then
                TW=Tr(i-1)
                TE=Tr(i+1)
                fT(i)=-aa*(TE-TP) + aa*(TP-TW) - ab*(rr-dr)*(TE+TP)**4 + ab*rr*(TP+TW)**4
            else
                TW=Tr(i-1)
                fT(n)= aa*(TP-TW) - ab*16*(rr-dr)*(TP)**4 + ab*rr*(TP+TW)**4
            end if
            rr=rr-dr
        end do
        fT=-fT
    end subroutine solve_fT

    subroutine solve_Jx(e, f, g, Tr)
        DOUBLE PRECISION, dimension(n) :: e, f, g, Tr
        DOUBLE PRECISION :: aa, ab, dr, rr, TW, TP, TE
        integer  :: i
    !!This subroutine creates and solve the Jacobian of our function
    !!Using Thomas TDMA algorithm

        !Variable initialization
        dr=r/n
        rr=r
        aa=1000/dr 
        ab=(11.34*10.0**(-5))/16
        do i = 1,n
            TP=Tr(i)
            if ( i==1 ) then
                TE=Tr(i+1)
                g(1)=-aa - 4*ab*(rr-dr)*(TE+TP)**3
                f(1)=3*aa - 4*ab*(rr-dr)*(TE+TP)**3
            else if ( i<n ) then
                TW=Tr(i-1)
                TE=Tr(i+1)
                g(i)=-aa - 4*ab*(rr-dr)*(TE+TP)**3
                f(i)=2*aa - 4*ab*(rr-dr)*(TE+TP)**3 + 4*ab*rr*(TP+TW)**3
                e(i)=-aa + 4*ab*rr*(TP+TW)**3 
            else
                TW=Tr(i-1)
                f(n)=aa - 64*ab*(rr-dr)*TP**3 + 4*ab*rr*(TP+TW)**3
                e(n)=-aa + 4*ab*rr*(TP+TW)**3 
            end if
            rr=rr-dr
        end do
    end subroutine solve_Jx
!*************************************************************************************
!*************************************************************************************
    subroutine solve_P(Tr, P1, P2)
        DOUBLE PRECISION, INTENT(OUT) :: P1, P2
        DOUBLE PRECISION, dimension(n) ::  qin, qout
        DOUBLE PRECISION, dimension(n) :: Tr
        DOUBLE PRECISION :: aa, ac, ab, TW, TP, rr, dr
        INTEGER :: i
    !!This subroutine solves the Power-in and Power-out, using our formulation

        !Variable initialization
        dr=r/n
        rr=r
        aa=0.001*pi
        ac=pi*11.34*10.0**(-8)   
        ab=aa*4*1000
        TW=Teg   !Temperature maintained at the edge
        do i = 1,n
            TP=Tr(i)
            if ( i==1 ) then
                qin(1)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(1)=ac*(rr**2-(rr-dr)**2)*TP**4
            else if ( i<n ) then
                TW=Tr(i-1)
                qin(i)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(i)=ac*(rr**2-(rr-dr)**2)*TP**4
            else
                TW=Tr(i-1)
                qin(n)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(n)=ac*(rr**2-(rr-dr)**2)*TP**4
            end if
            rr=rr-dr
        end do
        P1=sum(qin)
        P2=2*sum(qout)!Multiplied by 2, because radiation occurs from the top and bottom.
    end subroutine solve_P

!*****************************************************************************************************
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
    !*****************************************************************************************************
end program newt_Jac_TDMA