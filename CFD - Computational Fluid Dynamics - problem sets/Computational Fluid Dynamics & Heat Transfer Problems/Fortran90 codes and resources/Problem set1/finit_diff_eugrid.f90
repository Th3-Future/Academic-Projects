!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on May 20, 2021.
!This program solves a linear/non-linear equation using finite difference method
!and was written as a solution to AE8112 PS1 q3b
!****************************End Header*****************************************
program finit_diff_eugrid
    !Variable declaration
    implicit none
    real, dimension(:,:), ALLOCATABLE:: Ag
    real, DIMENSION(:), ALLOCATABLE:: bg, Ti, Te
    real, PARAMETER :: L=1.0 !Domain length parameter definition
    real  :: Alp, Linf
    integer :: N, Nd
    N=81   !The number of discretized points
    Alp=0.7
    Nd=19  !The number of times we are dividing Alpha (it must be an odd value)
    open(1, file = 'diffu_error.csv', status = 'unknown')
do while(Alp<=1.3 .and. Alp/=1)  !This will run for 0.7<α<1.3.
   call creat_uTDMA(Ag,bg,Ti,Te) !call the subroutine to solve for unequispaced grid
    Linf=maxval(abs(Ti(1:n)-Te(1:n)))
    !outputing the results
    print 3, 'For α =',Alp
    write(*,1)  "T(x) = ",Ti
    write(*,2) "L∞_u = ",Linf
    !writing one of the errors (Linf) as a function of Alpha (α) to excel file for ploting
    write(1,*) Alp, Linf
    Alp=Alp+(1.3-0.7)/Nd
end do
call creat_eTDMA(Ti,Te)    !call the subroutine to solve for equispaced grid
Linf=maxval(abs(Ti(1:n)-Te(1:n)))
write(*,2) "L∞_e= ", Linf
1 format(a6,81f9.3)
2 format(a7,9f9.5)
3 format(a10,f5.2)
close(1)

contains
!*****************************************************************************************************  
subroutine creat_uTDMA(Ag1,bg1,Ti1,Te1)
    implicit none
    real, dimension(:,:), ALLOCATABLE, INTENT(OUT):: Ag1
    real, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: bg1, Ti1, Te1
    real :: dx1
    real  :: aw, ae, ap, aa, xi, b
    integer :: i, Ns
    !This subroutine generates the tri-diagonal matric (TDMA) for unequal grid point

    ALLOCATE(Ag1(n,n),bg1(n),Ti1(n),Te1(n)) !array allocation
    !Initializing Variables
    Ag1=0 
    xi=0    
    Ns=n-1   !this is the number of segment
    dx1=L*(1-Alp)/(1-Alp**Ns)  
    aa=1-dx1**3  !This is also noticed to yield same result as, aa=1/dx1, (they are enterchangable)
    do i = 1, n  !This loop Matrix composition of the given problem
        aw=2/dx1
        ae=2/(Alp*dx1)
        ap=aw+ae
        b=100*(xi**2)*(dx1+Alp*dx1)
        !for the the first gridpoint, applying B.Cs T(0)=1
        if ( i==1 ) then
            Ag1(i,i)=-aa
            Ag1(i,i+1)=aa
            bg1(i)=0
        !for the the intermediate gridpoint
        else if ( i<n ) then
            Ag1(i,i-1)=aw
            Ag1(i,i)=-ap
            Ag1(i,i+1)=ae
            bg1(i)=-b
            dx1=Alp*dx1
        !for the the last gridpoint, applying B.Cs T(2pi)=1
        else 
            Ag1(n,n)=1
            bg1(n)=0
            dx1=Alp*dx1
        end if
        Te1(i)=(25.0/3)*(1-xi**4)
        xi=xi+dx1
    end do
    call tdma(Ag1,bg1,Ti1)
end subroutine creat_uTDMA  
!*****************************************************************************************************

!*****************************************************************************************************
subroutine creat_eTDMA(Ti1,Te1)
    implicit none
    real, dimension(:,:), ALLOCATABLE :: Ag1
    real, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: Ti1, Te1
    real, DIMENSION(:), ALLOCATABLE:: bg1
    real:: dx1
    real  :: a, ap, aa, xi, b
    integer :: i, Ns
    !This subroutine generates the tri-diagonal matric (TDMA) for equal grid spacing

    ALLOCATE(Ag1(n,n),bg1(n),Ti1(n),Te1(n)) !array allocation
    !Initializing Variables
    Ag1=0 
    xi=0    
    Ns=n-1   !this is the number of segment
    dx1=L/(n-1) 
    aa=1/dx1  !This is also noticed to yield same result as, aa=1/dx1, (they are enterchangable)
    do i = 1, n  !This loop Matrix composition of the given problem
        a=1/dx1**2
        ap=a+a
        b=100*(xi**2)
        !for the the first gridpoint, applying B.Cs T(0)=1
        if ( i==1 ) then
            Ag1(i,i)=-aa
            Ag1(i,i+1)=aa
            bg1(i)=0
        !for the the intermediate gridpoint
        else if ( i<n ) then
            Ag1(i,i-1)=a
            Ag1(i,i)=-ap
            Ag1(i,i+1)=a
            bg1(i)=-b
        !for the the last gridpoint, applying B.Cs T(2pi)=1
        else 
            Ag1(n,n)=1
            bg1(n)=0
        end if
        Te1(i)=(25.0/3)*(1-xi**4)
        xi=xi+dx1
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

end program finit_diff_eugrid