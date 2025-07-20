!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS1 q3b
!****************************End Header*****************************************
program Vfinit_pBi_CGSTAB
    !Variable declaration
    implicit none
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE:: Ag
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: bg, Ti, Te
    DOUBLE PRECISION, PARAMETER :: Pi = 4*atan(1.0) !pi parameter definition
    DOUBLE PRECISION  :: L1, L2, Linf, Rs, dx
    integer :: N, ii
    N=64 !Number of control volumes
    ALLOCATE(Ag(n,n),bg(n),Ti(n),Te(n)) !array allocation
   call creat_eTDMA(Ag,bg,Te,dx)
   ! for the initial guess
    do ii = 1, n
        Ti(ii)=bg(ii)/Ag(ii,ii)
    end do
    call Bi_CGSTAB_P(Ag,bg,Ti,ii)
   !Calculating the error terms
    L1=sum(dabs(Ti(1:n)-Te(1:n)))
    L2=sum((Ti(1:n)-Te(1:n))**2)
    Linf=maxval(dabs(Ti(1:n)-Te(1:n)))

    Rs=sqrt(sum((bg-matmul(Ag,Ti))**2))/N !Calculates the linear system residual term
    !Printing our results
    print 3, 'For number Iteration i =',ii
    write(*,1)  "T(x) = ",Ti
    write(*,2) "L1 = ", L1
    write(*,2) "L2 = ", L2
    write(*,2) "L∞ = ", Linf
    write(*,2) "R = ", Rs
    1 format(a6,64f8.3)
    2 format(a7,E10.3)
    3 format(a30,i3)

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
subroutine Bi_CGSTAB_P(A,b,x,i)
    implicit none
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: x
    DOUBLE PRECISION, dimension(n,n) :: A, K
    DOUBLE PRECISION, dimension(n), INTENT(IN):: b
    DOUBLE PRECISION, DIMENSION(n) ::   p, r, r0, y, z, v, s, t
    DOUBLE PRECISION :: rho0, rho, w, alpha, beta, r_check
    INTEGER :: i
!!This subroutine solves a tri-diagonal linear system using The Preconditioned Bi_CGSTAB Algorithm
!!by Van Der Vorst

    !Variable initialization
    r0=b-matmul(A,x); r=r0
    w=1; alpha=1; rho0=1; r_check=1
    v=0; p=0;  k=0
    !For the inverse of K
    do i = 1, n
        k(i,i)=1/A(i,i) 
    end do
    i=0
    open(1, file = 'PS2_Q3b.txt', status = 'unknown')
    !main algorithm loop
    do while (r_check>=0.663E-6)
        i=i+1
        rho=dot_product(r0,r)
        beta=(rho/rho0)*(alpha/w)
        rho0=dot_product(r0,r)
        p=r+beta*(p-w*v)
        y=matmul(K,p)
        v=matmul(A,y)
        alpha=rho/dot_product(r0,v)
        s=r-alpha*v
        z=matmul(K,s)
        t=matmul(A,z)
        w=dot_product(matmul(K,t),matmul(K,s))/dot_product(matmul(K,t),matmul(K,t))
        x=x+alpha*y+w*z
        r=s-w*t 
        r_check=norm2(r)
        write(1,*) i, r_check !writes the results to .txt file
    end do
    close(1)
end subroutine Bi_CGSTAB_P
!*************************************************************************************

end program Vfinit_pBi_CGSTAB