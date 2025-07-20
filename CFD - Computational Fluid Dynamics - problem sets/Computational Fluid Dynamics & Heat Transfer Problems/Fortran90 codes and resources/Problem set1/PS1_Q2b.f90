!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on May 20, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS1 q2b
!****************************End Header*****************************************
program finit_vol
    !Variable declaration
    implicit none
    real, dimension(:,:), ALLOCATABLE:: Ag
    real, DIMENSION(:), ALLOCATABLE:: bg, Ti, Te
    real, PARAMETER :: Pi = 4*atan(1.0) !pi parameter definition
    real  :: L1, L2, Linf, dx
    integer :: N
    N=8
    open(1, file = 'vol_error.csv', status = 'unknown') 
do while(n<=64)  !This will run for N=8,16,32 and 64
   call creat_eTDMA(Ag,bg,Ti,Te,dx)
   !Calculating the error terms
    L1=sum(abs(Ti(1:n)-Te(1:n)))
    L2=sum((Ti(1:n)-Te(1:n))**2)
    Linf=maxval(abs(Ti(1:n)-Te(1:n)))
    !outputing the results
    print 3, 'For number of grid points, n =',n
    write(*,1)  "T(x) = ",Ti
    write(*,2) "L1 = ", L1
    write(*,2) "L2 = ", L2
    write(*,2) "L∞ = ", Linf
    1 format(a6,64f8.3)
    2 format(a7,8f8.5)
    3 format(a30,i3)
    !writing one of the errors (L2) as a function of grid spacing (dx) to excel file for ploting
    write(1,*) dx, L2
    N=n+n
end do
close(1)

contains
!*****************************************************************************************************  
subroutine creat_eTDMA(Ag1,bg1,Ti1,Te1,dx1)
    implicit none
    real, dimension(:,:), ALLOCATABLE, INTENT(OUT):: Ag1
    real, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: bg1, Ti1, Te1
    real, INTENT(OUT) :: dx1
    real  :: aa, ap, xi
    integer :: i
    !This subroutine generates the tri-diagonal matric (TDMA) using our derived equations

    ALLOCATE(Ag1(n,n),bg1(n),Ti1(n),Te1(n)) !array allocation
    !Initializing Variables
    Ag1=0 
    xi=0
        dx1=2*pi/n
        aa=1/(dx1)
        ap=aa+aa
    do i = 1, n  !This loop Matrix composition of the given problem
        xi=xi+dx1
        !for the the first gridpoint, applying B.Cs T(0)=1
        if ( i==1 ) then
            Ag1(i,i)=-aa-2*aa
            Ag1(i,i+1)=aa
            bg1(i)=-2*aa-sin(xi)
            Te1(i)=cos(dx1/2) !This is T_exact from our analytical solution
        !for the the intermediate gridpoint
        else if ( i<n ) then
            Ag1(i,i-1)=aa
            Ag1(i,i)=-ap
            Ag1(i,i+1)=aa
            bg1(i)=sin(xi-dx1)-sin(xi)
            Te1(i)=cos(xi-dx1/2)
        !for the the last gridpoint, applying B.Cs T(2pi)=1
        else 
            Ag1(n,n-1)=aa
            Ag1(n,n)=-2*aa-aa
            bg1(n)=-2*aa+sin(xi-dx1)
            Te1(n)=cos(xi-dx1/2)
        end if
    end do
    call tdma(Ag1,bg1,Ti1)
end subroutine creat_eTDMA  
!*****************************************************************************************************

!*****************************************************************************************************
subroutine tdma(A,b1,x)
    implicit none
    real, dimension(n,n), INTENT(IN):: A
    real, dimension(n), INTENT(IN):: b1
    real, DIMENSION(n), INTENT(OUT) :: x
    real, DIMENSION(n):: b, e, f, g
    INTEGER :: k
    !This subroutine solves a tri-diagonal linear system using the Thomas Algorithm

    b=b1
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