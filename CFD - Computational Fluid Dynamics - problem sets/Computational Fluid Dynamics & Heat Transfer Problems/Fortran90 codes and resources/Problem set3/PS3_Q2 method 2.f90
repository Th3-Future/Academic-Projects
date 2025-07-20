!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves an unsteady linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS3 q2
!****************************End Header*****************************************
program unsteady_Vfinite
    implicit none
    !Variable declaration
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: Ti
    DOUBLE PRECISION, PARAMETER :: rh=1716, Cp=4817, lamda=14.6, ha=472, hs=36.4
    DOUBLE PRECISION, PARAMETER :: L=0.05, tk=0.0035 !plate dimension
    DOUBLE PRECISION :: dx, dt, tt, Tinf, tols, start, finish
    INTEGER :: N, ii, i1, j1, l1
    call cpu_time(start) !gets the start time, for timing purpose
    open(1, file = 'PS3_Q22.txt', status = 'unknown') 
    N=3 !Number of Control Volume along x&y direction
    dx=L/n
    dt=0.01     !given timestep
    Tinf=298.15 !Ambient temperature
    ALLOCATE(Ti(n*n))
    do ii = 1,4
        Ti=298.15 !initial temperature of the plate
        if (ii==1) then
            tols=5E-3
        elseif (ii==2) then
            tols=5E-5
        elseif (ii==3) then
            tols=5E-8
        else
            tols=5E-10
        end if
        call unsteady(Ti, tt) !solves the unsteady system
        !Printing results
        write(1,2) "Tolerance = ",tols, "t_total = ",tt-dt,"sec"
        j1=n*n-(n-1)
        l1=n*n
        do i1 = j1, 1, -n
            !write(*,1) Ti(i1:l1)
            write(1,1) Ti(i1:l1)
            l1=l1-n
        end do
        Print *, "t_total = ", tt-dt,"sec"
    end do
    call cpu_time(finish) !gets the end time
    write(*,3) "Program Execution Time = ", (finish-start)/60,"min"
    1 format(80f8.3)
    2 format(a15,E12.1,a15,E17.3,a4)
    3 format(a25,f7.3)
    close(1)

    contains
!*************************************************************************************
    subroutine unsteady(Tp, t)
        DOUBLE PRECISION, dimension(n*n) :: d, e, f, g, h, b, Tp, xi, T_old
        DOUBLE PRECISION :: aa, ap, bb, k1, k2, k3, alpha, t, b_simp, Tnorm
        integer  :: i, j, k
    !!This subroutine solves unsteady thermal problem using O'Brien's method

        
        !Variable initialization
        d=0; e=0; f=0; g=0; h=0
        t=0
        Tnorm=1
        alpha=lamda/(rh*Cp) !Thermal diffusivity
        aa=alpha*dt
        ap=(dx**2)+4*alpha*dt+(2*hs*dt*dx**2)/(rh*cp*tk)
        k1=ha*Tinf/((lamda/dx)+(ha/2))
        k2=((ha/2)-(lamda/dx))/((lamda/dx)+(ha/2))
        k3=846.3
        b_simp=(2*hs*dt*Tinf*dx**2)/(rh*cp*tk)
        do while (Tnorm>tols) !Loop for time step
            j=n+1
            k=n+n
            T_old=Tp
            do i = 1,n*n !Loop for control volume
                bb=Tp(i)*(dx**2)+b_simp
                !The 9 if-statements below are for the 9 regions of our CV formulation
                if ( i==1 ) then !bottom left region
                    f(i)=-(aa*k2+aa+ap)
                    g(i+1)=aa
                    h(i+n)=aa
                    b(i)=-(bb+aa*k3+aa*k1)
                    xi(i)=dx/2
                else if ( i < n ) then        !bottom region
                    e(i-1)=aa
                    f(i)=-(aa+ap)
                    g(i+1)=aa
                    h(i+n)=aa
                    b(i)=-(bb+aa*k3)
                    xi(i)=xi(i-1)+dx
                else if ( i==n ) then         !bottom right region
                    e(i-1)=aa
                    f(i)=-(aa*k2+aa+ap)
                    h(i+n)=aa
                    b(i)=-(bb+aa*k3+aa*k1)
                    xi(i)=xi(i-1)+dx
                else if ( i==n*n-(n-1) ) then !top left region
                    d(i-n)=aa
                    f(i)=-(2*aa*k2+ap)
                    g(i+1)=aa
                    b(i)=-(bb+2*aa*k1)
                else if ( i==n*n ) then       !top right region
                    d(i-n)=aa
                    e(i-1)=aa
                    f(i)=-(2*aa*k2+ap)
                    b(i)=-(bb+2*aa*k1)
                else if ( i==j ) then         !left region
                    d(i-n)=aa
                    f(i)=-(aa*k2+ap)
                    g(i+1)=aa
                    h(i+n)=aa
                    b(i)=-(bb+aa*k1)
                    j=j+n
                else if ( i==k ) then         !right region
                    d(i-n)=aa
                    e(i-1)=aa
                    f(i)=-(aa*k2+ap)
                    h(i+n)=aa
                    b(i)=-(bb+aa*k1)
                    k=k+n
                else if ( i < n*n-n ) then    !Interior region
                    d(i-n)=aa
                    e(i-1)=aa
                    f(i)=-ap
                    g(i+1)=aa
                    h(i+n)=aa
                    b(i)=-bb
                else                          !top region
                    d(i-n)=aa
                    e(i-1)=aa
                    f(i)=-(aa*k2+ap)
                    g(i+1)=aa
                    b(i)=-(bb+aa*k1)
                end if  
            end do
            !For computational efficiency, the A tri-diagonal matrix is split to 3 vectors
            call Bi_CGSTAB_P(d,e,f,g,h,b,Tp) !Solves the linear system at each time-step
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
    subroutine Bi_CGSTAB_P(d,e,f,g,h,b,x)
        implicit none
        DOUBLE PRECISION, DIMENSION(n*n), INTENT(OUT) :: x
        DOUBLE PRECISION, dimension(n*n), INTENT(IN):: d, e, f, g, h, b
        DOUBLE PRECISION, DIMENSION(n*n) ::   p, r, r0, y, z, v, s, t, k, ks, kt, wat
        DOUBLE PRECISION :: rho0, rho, w, alpha, beta, r_check, tol
        INTEGER :: i
    !!This subroutine solves a penta-diagonal linear system using The Preconditioned Bi_CGSTAB Algorithm
    !!by Van Der Vorst
    
        !Variable initialization
        wat=cshift(g(1:n*n)*x(1:n*n),1)
        r0=b-(cshift(d(1:n*n)*x(1:n*n),-n)+cshift(e(1:n*n)*x(1:n*n),-1)+f(1:n*n)*x(1:n*n)&
            +cshift(g(1:n*n)*x(1:n*n),1)+cshift(h(1:n*n)*x(1:n*n),n))
        r=r0
        w=1; alpha=1; rho0=1; r_check=1
        v=0; p=0;  k=0
        tol=10E-9
        !For the inverse of K
        do i = 1, n*n
            k(i)=1/f(i) 
        end do
        i=0
        !main algorithm loop
        do while (r_check>tol)
            i=i+1
            rho=dot_product(r0,r)
            beta=(rho/rho0)*(alpha/w)
            rho0=dot_product(r0,r)
            p=r+beta*(p-w*v)
            y=k(1:n*n)*p(1:n*n)
            v=(cshift(d(1:n*n)*y(1:n*n),-n)+cshift(e(1:n*n)*y(1:n*n),-1)+f(1:n*n)*y(1:n*n)&
                +cshift(g(1:n*n)*y(1:n*n),1)+cshift(h(1:n*n)*y(1:n*n),n))
            alpha=rho/dot_product(r0,v)
            s=r-alpha*v
            z=k(1:n*n)*s(1:n*n)
            t=(cshift(d(1:n*n)*z(1:n*n),-n)+cshift(e(1:n*n)*z(1:n*n),-1)+f(1:n*n)*z(1:n*n)&
                +cshift(g(1:n*n)*z(1:n*n),1)+cshift(h(1:n*n)*z(1:n*n),n))
            ks=(cshift(d(1:n*n)*s(1:n*n),-n)+cshift(e(1:n*n)*s(1:n*n),-1)+f(1:n*n)*s(1:n*n)&
                +cshift(g(1:n*n)*s(1:n*n),1)+cshift(h(1:n*n)*s(1:n*n),n))
            kt=(cshift(d(1:n*n)*t(1:n*n),-n)+cshift(e(1:n*n)*t(1:n*n),-1)+f(1:n*n)*t(1:n*n)&
                +cshift(g(1:n*n)*t(1:n*n),1)+cshift(h(1:n*n)*t(1:n*n),n))
            w=dot_product(kt,ks)/dot_product(kt,kt)
            x=x+alpha*y+w*z
            r=s-w*t 
            r_check=norm2(r)
        end do
    end subroutine Bi_CGSTAB_P
!*************************************************************************************
end program unsteady_Vfinite