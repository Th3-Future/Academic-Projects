!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves an unsteady linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS3 q2
!****************************End Header*****************************************
program unsteady_Vfinite
    implicit none
    !Variable declaration
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: Ti
    DOUBLE PRECISION, PARAMETER :: rh=1716, Cp=4817, lamda=14.6, h=472, hs=36.4
    DOUBLE PRECISION, PARAMETER :: L=0.05, tk=0.0035 !plate dimension
    DOUBLE PRECISION :: dx, dt, tt, Tinf
    INTEGER :: N, i1, j1, l1
    open(1, file = 'PS3_Q2.txt', status = 'unknown') 
    N=10 !Number of Control Volume along x&y direction
    dx=L/n
    dt=0.01     !given timestep
    ALLOCATE(Ti(n*n))
    Tinf=298.15 !Ambient temperature
    Ti=298.15 !initial temperature of the plate
    call unsteady(Ti, tt) !solves the unsteady system
    !Printing results
    j1=n*n-(n-1)
    l1=n*n
    do i1 = j1, 1, -n
        write(*,1) Ti(i1:l1)
        write(1,1) Ti(i1:l1)
        l1=l1-n
    end do
    Print *, "t_total = ", tt-dt,"sec"

    1 format(80f8.3)
    close(1)
    contains
!*************************************************************************************
    subroutine unsteady(Tp, t)
        DOUBLE PRECISION, dimension(n*n,n*n) :: A
        DOUBLE PRECISION, dimension(n*n) :: b, Tp, xi, T_old
        DOUBLE PRECISION :: aa, ap, bb, k1, k2, k3, alpha, t, Tnorm, tol
        integer  :: i, j, k
    !!This subroutine solves unsteady thermal problem using O'Brien's method

        
        !Variable initialization
        A=0
        t=0
        tol=5E-3
        Tnorm=1
        alpha=lamda/(rh*Cp) !Thermal diffusivity
        aa=alpha*dt
        ap=(dx**2)+4*alpha*dt+(2*hs*dt*dx**2)/(rh*cp*tk)
        k1=h*Tinf/((lamda/dx)+(h/2))
        k2=((h/2)-(lamda/dx))/((lamda/dx)+(h/2))
        k3=846.3
        do while (Tnorm>tol) !Loop for time step
            j=n+1
            k=n+n
            do i = 1,n*n !Loop for control volume
                bb=Tp(i)*(dx**2)+(2*hs*dt*Tinf*dx**2)/(rh*cp*tk)
                T_old=Tp
                !The 9 if-statements below are for the 9 regions of our CV formulation
                if ( i==1 ) then !bottom left region
                    A(i,i)=-(aa*k2+aa+ap)
                    A(i,i+1)=aa
                    A(i,i+n)=aa
                    b(i)=-(bb+aa*k3+aa*k1)
                    xi(i)=dx/2
                else if ( i < n ) then        !bottom region
                    A(i,i-1)=aa
                    A(i,i)=-(aa+ap)
                    A(i,i+1)=aa
                    A(i,i+n)=aa
                    b(i)=-(bb+aa*k3)
                    xi(i)=xi(i-1)+dx
                else if ( i==n ) then         !bottom right region
                    A(i,i-1)=aa
                    A(i,i)=-(aa*k2+aa+ap)
                    A(i,i+n)=aa
                    b(i)=-(bb+aa*k3+aa*k1)
                    xi(i)=xi(i-1)+dx
                else if ( i==n*n-(n-1) ) then !top left region
                    A(i,i-n)=aa
                    A(i,i)=-(2*aa*k2+ap)
                    A(i,i+1)=aa
                    b(i)=-(bb+2*aa*k1)
                else if ( i==n*n ) then       !top right region
                    A(i,i-1)=aa
                    A(i,i-n)=aa
                    A(i,i)=-(2*aa*k2+ap)
                    b(i)=-(bb+2*aa*k1)
                else if ( i==j ) then         !left region
                    A(i,i-n)=aa
                    A(i,i)=-(aa*k2+ap)
                    A(i,i+1)=aa
                    A(i,i+n)=aa
                    b(i)=-(bb+aa*k1)
                    j=j+n
                else if ( i==k ) then         !right region
                    A(i,i-1)=aa
                    A(i,i-n)=aa
                    A(i,i)=-(aa*k2+ap)
                    A(i,i+n)=aa
                    b(i)=-(bb+aa*k1)
                    k=k+n
                else if ( i < n*n-n ) then    !Interior region
                    A(i,i-1)=aa
                    A(i,i-n)=aa
                    A(i,i)=-ap
                    A(i,i+1)=aa
                    A(i,i+n)=aa
                    b(i)=-bb
                else                          !top region
                    A(i,i-1)=aa
                    A(i,i-n)=aa
                    A(i,i)=-(aa*k2+ap)
                    A(i,i+1)=aa
                    b(i)=-(bb+aa*k1)
                end if  
            end do
            !For computational efficiency, the A tri-diagonal matrix is split to 3 vectors
            call Bi_CGSTAB_P(A,b,Tp) !Solves the linear system at each time-step
            ! if ( t==0 ) then
            !     write(1,3) "Time(sec)", "T(0,t)", (xi(j), j=1,n) !writes the result's title to .txt file
            ! end if
            !write(1,4) t*1D-5, Tww, (Tp(j), j=1,n) !writes the results to .txt file
            Tnorm=(sqrt(sum((T_old(1:n*n)-Tp(1:n*n))**2)))/n*n
            t=t+dt
        end do
        !3 format(a9,3x,a6,250(2x,f6.3))
        !4 format(f9.7,3x,f7.3,250(1x,f7.3))
    end subroutine unsteady
!*************************************************************************************
!*************************************************************************************
    subroutine Bi_CGSTAB_P(A,b,x)
        implicit none
        DOUBLE PRECISION, DIMENSION(n*n), INTENT(OUT) :: x
        DOUBLE PRECISION, dimension(n*n,n*n) :: A, K
        DOUBLE PRECISION, dimension(n*n), INTENT(IN):: b
        DOUBLE PRECISION, DIMENSION(n*n) ::   p, r, r0, y, z, v, s, t
        DOUBLE PRECISION :: rho0, rho, w, alpha, beta, r_check, tol
        INTEGER :: i
    !!This subroutine solves a penta-diagonal linear system using The Preconditioned Bi_CGSTAB Algorithm
    !!by Van Der Vorst
    
        !Variable initialization
        r0=b-matmul(A,x); r=r0
        w=1; alpha=1; rho0=1; r_check=1
        v=0; p=0;  k=0
        tol=10E-9
        !For the inverse of K
        do i = 1, n*n
            k(i,i)=1/A(i,i) 
        end do
        i=0
        !main algorithm loop
        do while (r_check>tol)
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
        end do
    end subroutine Bi_CGSTAB_P
!*************************************************************************************
end program unsteady_Vfinite