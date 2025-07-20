!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 20, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS3 q1a
!****************************End Header*****************************************
program finit_vol
    !Variable declaration
    implicit none
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE:: Ag
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: bg, Phi_x, Phi_ex, x1
    DOUBLE PRECISION, PARAMETER :: L=1, Pe=1.09639004471!length of the domain and Peclet number
    DOUBLE PRECISION  :: L2, dx1, alpha
    integer :: N, ii, jj
    open(1, file = 'PS3_Q1a_error.txt', status = 'unknown') 
    open(2, file = 'PS3_Q1a_soln.txt', status = 'unknown')
    open(3, file = 'PS3_Q1b.txt', status = 'unknown')
    open(4, file = 'PS3_Q1b_pe.txt', status = 'unknown')
    write(1,5) 'L2'
    write(3,7) 'L2'
    write(4,5) 'L2'
    do ii = 1, 3
        N=32
        if ( ii==1 ) then !For central difference scheme
            alpha=0
        elseif (ii==2) then !For upwinded scheme
            alpha=1
        else 
            alpha=pe**2/(pe**2+5)
        end if
        do while(n<=128)  !This will run for N=32, 64 and 128
            ALLOCATE(Ag(n,n),bg(n),Phi_x(n),Phi_ex(n),x1(n))   
            call creat_eTDMA(Ag,bg,Phi_x,Phi_ex,dx1,x1)
            !Calculating L2 error term
            L2=sum((Phi_x(1:n)-Phi_ex(1:n))**2)/N
            !outputing the results
            print 3, 'For number of grid points, n =',n
            print 4, 'For α =',alpha
            write(*,1)  "φ(x) = ",Phi_x
            write(*,1)  "φ_exact = ",Phi_ex
            write(*,2) "L2 = ", L2
            1 format(a6,128f8.3)
            2 format(a7,E9.3)
            3 format(a30,i7)
            4 format(a7,f7.5)
            if ( n==64 ) then
                if ( alpha==0 ) then
                    write(2,*) 'solution for α =', alpha
                elseif (alpha==1) then
                    write(2,*) 'solution for α =', alpha
                else
                    write(2,*) 'solution for α =', alpha
                end if
                do jj = 1, n
                    write(2,*) x1(jj), Phi_x(jj)
                end do
            end if
            !writing the nine L2 error terms
            write(1,6) 'For No. of CVs. = ',N,', And α = ',alpha, L2
            write(3,8) 'Δx = ',dx1,', And α = ',alpha, L2
            if ( alpha==0 .or. alpha==1 ) then
                write(4,10) 'Δx = ',log(dx1),', And α = ',alpha, log(L2)
            end if
            N=n+n
            DEALLOCATE(Ag, bg, Phi_x, Phi_ex, x1)
        end do
    end do
    5 format(40x,a2)
    7 format(30x,a2)
    6 format(a19,i3,a11,f5.3,2x,E9.3)
    8 format(a6,f7.5,a11,f5.3,2x,E9.3)
    10 format(a6,E12.5,a11,f5.3,2x,E10.3)
    close(1);close(2);close(3);close(4)

contains
!*****************************************************************************************************  
subroutine creat_eTDMA(A,b,phi,phi_e,dx,x)
    implicit none
    DOUBLE PRECISION, dimension(n,n), INTENT(OUT):: A
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT):: b, phi, phi_e, x
    DOUBLE PRECISION, INTENT(OUT) :: dx
    DOUBLE PRECISION  :: ae, ap, aw, xi, af, al, ab
    integer :: i
    !This subroutine generates the tri-diagonal matric using our derived equations

    !Initializing Variables
    A=0; b=0
    dx=L/n
    xi=dx/2
    ab=((2/dx)+(Pe/L))
    af=(3/dx)+(Pe/(2*L))*(1+alpha)
    al=(3/dx)-(Pe/(2*L))*(1-alpha)
    ae=(1/dx)-(Pe/(2*L))*(1-alpha)
    ap=(2/dx)+(Pe*alpha/L)
    aw=(1/dx)+(Pe/(2*L))*(1+alpha)
    do i = 1, n  !This loop Matrix composition of the given problem 
        !for the the first gridpoint, applying B.Cs φ(0)=1
        if ( i==1 ) then
            A(i,i)=-(aw+ap)
            A(i,i+1)=ae
            b(i)=-(2*aw)
        !for the the intermediate gridpoint
        else if ( i<n ) then
            A(i,i-1)=aw
            A(i,i)=-ap
            A(i,i+1)=ae
            b(i)=0
        !for the the last gridpoint, applying B.Cs φ(L)=0
        else 
            A(n,n-1)=aw
            A(n,n)=-(ae+ap)
            b(n)=0
        end if
        phi_e(i)=(exp(pe)-exp(pe*xi/L))/(exp(pe)-1) !This is φ_exact from our analytical solution
        x(i)=xi
        xi=xi+dx
    end do
    call tdma(A,b,phi)
end subroutine creat_eTDMA  
!*****************************************************************************************************

!*****************************************************************************************************
subroutine tdma(A,b1,x)
    implicit none
    DOUBLE PRECISION, dimension(n,n), INTENT(IN):: A
    DOUBLE PRECISION, dimension(n), INTENT(IN):: b1
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: x
    DOUBLE PRECISION, DIMENSION(n):: b, e, f, g
    INTEGER :: k
    !This subroutine solves a tri-diagonal linear system using the Thomas Algorithm

    b=b1
    x=0
   !extracting e, f and g array
    do k=1,n-1
        e(k+1)=A(k+1,k)
        g(k)=A(k,k+1)
    end do 
    do  k = 1,n   
        f(k)=A(k,k)
    end do
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

end program finit_vol