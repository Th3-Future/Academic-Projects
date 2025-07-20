!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS1 q3a
!****************************End Header*****************************************
program Vfinite_Gauss_seid
    !Variable declaration
    implicit none
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE:: Ag
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: bg, Ti, Te
    DOUBLE PRECISION, PARAMETER :: Pi = 4*atan(1.0) !pi parameter definition
    DOUBLE PRECISION  :: L1, L2, Linf, R, dx
    integer :: N, ii, m
    N=8 !First number of control volumes
    open(1, file = 'PS2_Q3a.txt', status = 'unknown')
do while(n<=64)  !This will run for N=8,16,32 and 64
    ALLOCATE(Ag(n,n),bg(n),Ti(n),Te(n)) !array allocation
   call creat_eTDMA(Ag,bg,Te,dx) !creates Tri-diagonal matric for equispaced grid
   ! for the initial guess
    do ii = 1, n
        Ti(ii)=bg(ii)/Ag(ii,ii)
    end do
   !For Gauss Siedel iteration
   do m = 1, 3000
        call gauss_siedel(Ag,bg,Ti)
        R=sqrt(sum((bg-matmul(Ag,Ti))**2))/N
        write(1,*) m, R
   end do
   
   close(1)
   !Calculating the error terms
    L1=sum(dabs(Ti(1:n)-Te(1:n)))
    L2=sum((Ti(1:n)-Te(1:n))**2)
    Linf=maxval(dabs(Ti(1:n)-Te(1:n)))
    !Calculating the linear system residual term
    
    print 3, 'For number of grid points, n =',n
    write(*,1)  "T(x) = ",Ti
    write(*,2) "L1 = ", L1
    write(*,2) "L2 = ", L2
    write(*,2) "L∞ = ", Linf
    write(*,2) "R = ", R
    1 format(a6,64f8.3)
    2 format(a7,E10.3)
    3 format(a30,i3)
    !writing one of the errors (L2) as a function of grid spacing (dx) to a .txt file
    write(1,*) dx, L2
    if(n<=64) deallocate(Ag,bg,Ti,Te) !deallocates the arrays for every new N
    N=n+n
end do

contains
!*************************************************************************************
subroutine creat_eTDMA(Ag1,bg1,Te1,dx1)
    implicit none
    DOUBLE PRECISION, dimension(n,n), INTENT(OUT):: Ag1
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT):: bg1, Te1
    DOUBLE PRECISION, INTENT(OUT) :: dx1
    DOUBLE PRECISION  :: aa, ap, xi
    integer :: i
    !!This subroutine generates the tri-diagonal matric using our derived equations

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
            bg1(i)=-2*aa-cos(dx1/2)*dx1
            Te1(i)=cos(dx1/2) !This is T_exact from our analytical solution
        !for the the intermediate gridpoint
        else if ( i<n ) then
            Ag1(i,i-1)=aa
            Ag1(i,i)=-ap
            Ag1(i,i+1)=aa
            bg1(i)=-cos(xi-dx1/2)*dx1
            Te1(i)=cos(xi-dx1/2)
        !for the the last gridpoint, applying B.Cs T(2pi)=1
        else 
            Ag1(n,n-1)=aa
            Ag1(n,n)=-2*aa-aa
            bg1(n)=-2*aa-cos(xi-dx1/2)*dx1
            Te1(n)=cos(xi-dx1/2)
        end if
    end do
    
end subroutine creat_eTDMA  
!*************************************************************************************
!*************************************************************************************
subroutine gauss_siedel(A,b,x)
    implicit none
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: x
    DOUBLE PRECISION, dimension(n,n), INTENT(IN):: A
    DOUBLE PRECISION, dimension(n), INTENT(IN):: b
    DOUBLE PRECISION, DIMENSION(n) :: x_m
    INTEGER :: i
!!This subroutine solves a tri-diagonal linear system using the Gauss Siedel Method
    do i = 1,N
        x_m(i)=(1/a(i,i))*(b(i)-sum(a(i,1:(i-1))*x_m(1:(i-1)))-sum(a(i,(i+1):N)*x((i+1):N)))
    end do
    x=x_m
end subroutine gauss_siedel
!*************************************************************************************

end program Vfinite_Gauss_seid